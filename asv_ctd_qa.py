#!/usr/bin/env python
# coding=utf-8
import argparse
import calendar
import csv
import datetime
import json
import logging
import os

import netCDF4 as nc
import numpy as np
import pandas as pd
from cf_units import Unit  # noqa
from compliance_checker.runner import CheckSuite, ComplianceChecker

import ncdump
import qc_plots
from ioos_qc import qartod
from ioos_qc.config import QcConfig
from ioos_qc.qartod import QartodFlags
from ioos_qc.stores import PandasStore
from ioos_qc.streams import PandasStream

# Keyword to look for indicating the start of the data column row
INIT_COLUMN_ROW = "Date / Time"

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(asctime)s : %(msecs)04d : %(name)s : %(levelname)s : %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class UnicodeReader:
    """A CSV reader which will iterate over lines in the CSV file "f",.

    which is encoded in the given encoding.
    """

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds) -> None:
        self.reader = csv.reader(f, dialect=dialect, **kwds)
        self.rows_read = 0

    def next(self) -> list:
        row = next(self.reader)
        return row

    def __iter__(self):
        return iter(self)


class FileParser:
    """Parse input file attributes."""

    def __init__(self, fname: str, header_rows: int, column_delimiter: str) -> None:
        self.filename = fname
        self.delimiter = column_delimiter
        self.headrows = header_rows
        self.datarows_start = self.headrows + 1
        self.cols = []
        self.filerows = 0
        self.datarows = 0
        self.encoding = None
        self.reader = self.file_reader()
        self.eval_file()

    def file_reader(self) -> UnicodeReader:
        """Determine file dialect and encoding, then create.

        and return a reader object
        """
        dialect = "utf-8"
        try:
            dialect = csv.Sniffer().sniff(open(self.filename).readline())
        except Exception:
            raise
        try:
            if self.encoding:
                reader = UnicodeReader(
                    open(self.filename),
                    dialect,
                    self.encoding,
                )
            else:
                reader = csv.reader(open(self.filename), dialect)
        except Exception:
            raise
        if self.headrows > 0:
            try:
                for row in range(self.headrows):
                    next(reader)
            except Exception:
                raise
        self.cols = next(reader)
        return reader

    def eval_file(self):
        """Determine number of total rows and data rows in file."""
        self.filerows = self.headrows + 1
        for datarow in self.reader:
            self.filerows += 1
            if len(datarow) > 0:
                self.datarows += 1
        return self

    def filestats(self) -> str:
        """Returns multi-line string of file information."""
        fstats = (
            "Filename: %s\nFile rows: %d\nHeader rows: %d\nColumn row: %d\nData rows: %d\n"
            % (
                self.filename,
                self.filerows,
                self.headrows,
                self.datarows_start,
                self.datarows,
            )
        )
        fstats += "Columns:\n"
        for col in self.cols:
            fstats += "  %s\n" % col
        return fstats


class Parameters(FileParser):
    def __init__(self, config, fname, header_rows, column_delimiter) -> None:
        super().__init__(fname, header_rows, column_delimiter)
        self.config = config
        self.reader = self.file_reader()
        self.data_array = []

    def parameter_dataframe(self) -> None:
        for (idx, row) in enumerate(self.reader):
            self.data_array.append(list(filter(None, row)))
        try:
            self.df = pd.DataFrame(self.data_array, columns=self.cols)
        except Exception:
            raise

        keys = {
            key: self.config["parameters"][key]["standard_name"]
            for key in self.config["parameters"]
        }
        logger.info("CREATING CF COMPLIANT PARAMETER NAMES")
        for key, value in keys.items():
            logger.info(f"  {key} -> {value}")
        self.df.rename(columns=keys, inplace=True)
        self.df["syntax_test"] = self.syntax_test()
        for key in keys:
            self.df[keys[key]] = self.df[keys[key]].astype(
                self.config["parameters"][key]["dtype"],
            )

    def syntax_test(self) -> None:
        # Start with everything as passing (1)
        flag_arr = np.ma.ones(self.datarows, dtype="uint8")
        for idx, row in enumerate(self.file_reader()):
            row = tuple(filter(None, row))
            if len(row) != len(self.cols):
                flag_arr[idx] = QartodFlags.FAIL
        if np.all(flag_arr == 4):
            raise Exception("All syntax tests have failed.")
        return flag_arr

    def convert_units(self):
        for variable in self.config["parameters"]:
            p = self.config["parameters"][variable]
            if "original_units" in p:
                if p["original_units"] != p["units"]:
                    v = p["standard_name"]
                    _from = p["original_units"]
                    _to = p["units"]
                    logger.info("CONVERTING {} UNITS {} --> {}".format(v, _from, _to))
                    try:
                        self.df[v] = Unit(_from).convert(
                            self.df[v].to_numpy(),
                            Unit(_to),
                        )
                    except Exception:
                        raise


class NetCDF:
    """Create a NetCDF object and provide methods for adding self-describing attributes."""

    def __init__(self, config, outfile: str, data_model="NETCDF4") -> None:
        self.filename = outfile
        self.data_model = data_model
        self.rootgrp = nc.Dataset(self.filename, "w", format=self.data_model)
        self.config = config
        self.metadata()
        self.create_dimensions()

    def metadata(self) -> None:
        for attr in self.config["attributes"]:
            if attr == "history":
                setattr(
                    self.rootgrp,
                    attr,
                    f"{self.config['attributes'][attr]} {datetime.datetime.now().isoformat()}",
                )
            else:
                setattr(self.rootgrp, attr, self.config["attributes"][attr])

    def create_dimensions(self) -> None:
        self.rootgrp.createDimension("time", None)

    def create_ancillary_variables(
        self,
        data: pd.DataFrame,
        parameter: str,
        variables: dict,
    ) -> None:
        # if df['ts'] = df.datetime.astype('int64')
        if variables["standard_name"] == "time":
            # dtype = "int64"
            dtype = "f8"
        else:
            dtype = variables["dtype"]
        v = self.rootgrp.createVariable(
            variables["standard_name"],
            np.dtype(dtype).char,
            ("time",),
        )
        v.original_name = parameter
        for attr in variables:
            if attr != "dtype":
                setattr(v, attr, variables[attr])
        if variables["standard_name"] == "time":
            v[:] = data.to_numpy().astype("datetime64[ns]")
        else:
            v[:] = data.to_numpy()

    def create_test_variables(
        self,
        data: pd.DataFrame,
        ancillary_variable: str,
        test_name: str,
        variables: dict,
    ) -> None:

        v = self.rootgrp.createVariable(
            test_name,
            np.int8,
            ("time",),
            fill_value=np.int8(9),
        )
        v.long_name = test_name
        v.standard_name = f"{ancillary_variable} status_flag"
        v.flag_values = np.array([1, 2, 3, 4, 9], dtype=np.int8)
        v.flag_meanings = "GOOD NOT_EVALUATED SUSPECT BAD MISSING"
        for attr in variables:
            if attr == "tinp" or variables[attr] is None:
                continue
            elif type(variables[attr]) == pd.Series:
                setattr(v, attr, variables[attr.to_numpy()])
            elif attr == "config":
                for member in variables[attr].members:
                    setattr(v, "config_fspan_minv", float(member.fspan.minv))
                    setattr(v, "config_fspan_maxv", float(member.fspan.maxv))
                    setattr(
                        v,
                        "config_tspan_minv",
                        float(
                            (member.tspan.minv - pd.Timestamp("1970-01-01"))
                            // pd.Timedelta("1s"),
                        ),
                    )
                    setattr(
                        v,
                        "config_tspan_maxv",
                        float(
                            (member.tspan.maxv - pd.Timestamp("1970-01-01"))
                            // pd.Timedelta("1s"),
                        ),
                    )
                    setattr(v, "config_vspan_minv", float(member.vspan.minv))
                    setattr(v, "config_vspan_maxv", float(member.vspan.maxv))
                    setattr(v, "config_zspan_minv", float(member.zspan.minv))
                    setattr(v, "config_zspan_maxv", float(member.zspan.maxv))
            else:
                setattr(v, attr, variables[attr])
        v[:] = data[test_name].to_numpy()


def create_dir(dir: str) -> None:
    """Execute the makedirs method on a given path.

    Arguments:
        :path: OS file path
    """
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except OSError:
            raise


def clparser() -> argparse.ArgumentParser:
    """Create a parser to handle input arguments and displaying.

    a script specific help message.
    """
    desc_msg = """Evaulate a data file of ASV CTD readings
        and apply quality assurence
        checks following QARTOD methods and assigning data
        quality flags as appropriate. Transform
        results into NetCDF format following IC standards."""
    parser = argparse.ArgumentParser(description=desc_msg)
    parser.add_argument(
        "config",
        help="Configuration JSON file.",
    )
    parser.add_argument(
        "input_file",
        help="Path to the input sensor data file.",
    )
    parser.add_argument(
        "output_dir",
        help="Directory for output files.",
    )
    parser.add_argument(
        "-l",
        "--log",
        action="store",
        dest="log",
        type=str,
        nargs="?",
        help="Path to a log file for script level logging",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        dest="plot",
        default=False,
        help="Create an HTML file containing plots of QC flags. Files are stored under a subdirectory of the specified output_dir.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="Control the amount of information to display.",
    )
    parser.add_argument(
        "-c",
        "--compliance",
        action="store_true",
        dest="compliance_check",
        default=False,
        help="Run IOOS compliance checker on compiled NetCDF file.",
    )
    return parser


def preamble(**kwargs) -> None:
    """Write a preamble to a log file consisting of key:value pairs."""
    logger.info("=" * 66)
    logger.info(
        "start time: %s %s"
        % (
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"),
            datetime.datetime.now().astimezone().tzinfo,
        ),
    )
    if len(kwargs) > 0:
        for arg in kwargs:
            logger.info(f"{arg}: {kwargs[arg]}")
    logger.info("=" * 66)


def current_time_window() -> tuple:
    """Calculate the current season time window.

    Returns a tuple with the following format:
        (season, start_time, end_time)
    """
    seasons = {
        "spring": range(3, 5),
        "summer": range(6, 8),
        "fall": range(9, 11),
    }
    today = datetime.datetime.now()
    current_season = None

    for season in seasons:
        if seasons[season].start <= today.month <= seasons[season].stop:
            current_season = season

    if not current_season:
        current_season = "winter"
        if today.month == 12:
            window_start = datetime.datetime.strptime(
                f"12 {today.year}",
                "%m %Y",
            )
            eom = calendar.monthrange(today.year + 1, 2)[1]
            window_end = datetime.datetime.strptime(
                f"2 {today.year + 1} 23:59:59",
                "%m %Y %H:%M:%S",
            )
        else:
            window_start = datetime.datetime.strptime(
                f"12 {today.year - 1}",
                "%m %Y",
            )
            eom = calendar.monthrange(today.year, 2)[1]
            window_end = datetime.datetime.strptime(
                f"2 {today.year} 23:59:59",
                "%m %Y %H:%M:%S",
            )
    else:
        window_start = datetime.datetime.strptime(
            f"{seasons[current_season].start} {today.year}",
            "%m %Y",
        )
        eom = calendar.monthrange(today.year, seasons[current_season].stop)[1]
        window_end = datetime.datetime.strptime(
            f"{eom} {seasons[current_season].stop} {today.year} 23:59:59",
            "%d %m %Y %H:%M:%S",
        )

    return (current_season, window_start, window_end)


def run_qc(
    config: dict,
    parameters: Parameters,
    ncfile: NetCDF,
) -> None:

    ioos_qc_config = config["ioos_qc"]
    for variable in ioos_qc_config["variables"]:
        variable_config = ioos_qc_config["variables"][variable]
        if "valid_range_test" in variable_config["axds"]:
            variable_config["axds"]["valid_range_test"]["valid_range"] = tuple(
                np.array(variable_config["axds"]["valid_range_test"]["valid_span"]),
            )
        if "density_inversion_test" in variable_config["qartod"]:
            variable_config["qartod"]["density_inversion_test"]["zinp"] = parameters.df[
                "depth"
            ].to_numpy()
        if "climatology_test" in variable_config["qartod"]:
            variable_config["qartod"]["climatology_test"]["inp"] = parameters.df[
                variable
            ].to_numpy()
            variable_config["qartod"]["climatology_test"]["tinp"] = parameters.df[
                "time"
            ]
            variable_config["qartod"]["climatology_test"]["zinp"] = parameters.df[
                "depth"
            ].to_numpy()
            season, window_start, window_end = current_time_window()

            cc = qartod.ClimatologyConfig()
            cc.add(
                tspan=tuple([window_start, window_end]),
                fspan=tuple(
                    np.array(
                        variable_config["qartod"]["climatology_test"]["fail_span"],
                    ),
                ),
                vspan=tuple(
                    np.array(
                        variable_config["qartod"]["climatology_test"]["suspect_span"],
                    ),
                ),
                zspan=tuple(
                    np.array(
                        variable_config["qartod"]["climatology_test"]["zspan"],
                    ),
                ),
            )
            variable_config["qartod"]["climatology_test"]["config"] = cc
        # Add global tests
        variable_config["qartod"].update(config["ioos_qc"]["global"]["qartod"])

    ps = PandasStream(parameters.df)
    qc = QcConfig(ioos_qc_config["variables"])
    qc_results = ps.run(qc)
    store = PandasStore(qc_results)
    results_store = store.save(write_data=False, write_axes=False)
    for variable in ioos_qc_config["variables"]:
        for test_type in ioos_qc_config["variables"][variable]:
            for test in ioos_qc_config["variables"][variable][test_type]:
                test_name = f"{variable}_{test_type}_{test}"
                logger.info(f"RUNNING TEST {test_name}")
                if ioos_qc_config["variables"][variable][test_type][test] is not None:
                    ncfile.create_test_variables(
                        results_store,
                        variable,
                        test_name,
                        ioos_qc_config["variables"][variable][test_type][test],
                    )
        # Apply primary QC flag which is an aggregate of all other QC tests
        flags = results_store[[c for c in results_store.columns if variable in c]]
        rollup_var = f"{variable}_rollup_qc"
        results_store[rollup_var] = qartod.qartod_compare(
            [flags[c].to_numpy() for c in flags],
        )
        ncfile.create_test_variables(
            results_store,
            variable,
            rollup_var,
            {},
        )


def get_num_headerrows(infile: str):
    """Open the input file and check for header rows.

    Searches the start of each line for a substring matching INIT_COLUMN_ROW
    """
    try:
        with open(infile) as f:
            lines = f.readlines()
    except OSError:
        raise
    header_row = None
    for idx, line in enumerate(lines):
        if line[: len(INIT_COLUMN_ROW)] == INIT_COLUMN_ROW:
            header_row = idx
            break
    if header_row is not None:
        return header_row
    else:
        raise Exception(
            "Could not determine start of data stream. Searched %i lines for string starting with %s."
            % (idx, INIT_COLUMN_ROW),
        )


def asv_ctd_qa(
    config: str,
    input_file: str,
    output_dir: str,
    plot: bool = False,
    verbose: bool = False,
    log=None,
    compliance_check: bool = False,
) -> None:
    """Run IOOS QC tests."""

    # Validate input arguments
    if not os.path.exists(config):
        raise OSError("The input configuration file does not exist: %s", config)
    if not os.path.splitext(config)[1] == ".json":
        raise OSError("The input configuration file is invalid: %s" % config)
    if not os.path.exists(input_file):
        raise OSError("The input data file does not exist: %s" % input_file)

    num_headerrows = get_num_headerrows(input_file)

    input_file_basefn = os.path.basename(input_file)
    create_dir(output_dir)
    out_ncfile = os.path.join(output_dir, f"{input_file_basefn}.nc")

    preamble(
        config=config,
        input_file=input_file,
        input_file_basefn=input_file_basefn,
        num_headerrows=num_headerrows,
        output_dir=output_dir,
        out_ncfile=out_ncfile,
        plot=plot,
        verbose=verbose,
    )

    try:
        with open(config, "r", encoding="utf-8") as f:
            config_json = json.load(f)
    except OSError:
        raise

    logger.info("EVAULATING INPUT FILE")
    parameters = Parameters(
        config=config_json,
        fname=input_file,
        header_rows=num_headerrows,
        column_delimiter="\t",
    )
    logger.info(parameters.filestats())
    logger.info("~" * 66)
    parameters.parameter_dataframe()

    parameters.convert_units()

    ncfile = NetCDF(config_json, out_ncfile)
    for parameter in config_json["parameters"]:
        logger.info(
            f"CREATING ATTRIBUTES FOR {config_json['parameters'][parameter]['standard_name']}",
        )
        ncfile.create_ancillary_variables(
            parameters.df[config_json["parameters"][parameter]["standard_name"]],
            parameter,
            config_json["parameters"][parameter],
        )

    run_qc(config_json, parameters, ncfile)

    # Pretty print NetCDF file information, also store sections as variables
    nc_attrs, nc_dims, nc_vars = ncdump.ncdump(
        ncfile,
        verb=True if verbose else False,
        log=log,
    )

    if compliance_check:
        logger.info("=" * 66)
        logger.info("RUNNING CF COMPLIANCE CHECKS")
        check_suite = CheckSuite()
        check_suite.load_all_available_checkers()
        compliance_dir = os.path.join(output_dir, "compliance_checks")
        convention = json.dumps(config_json["attributes"]["Conventions"]).strip('"')
        create_dir(compliance_dir)
        ComplianceChecker.run_checker(
            ds_loc=out_ncfile,
            checker_names=[convention.replace("-", ":").lower()],
            verbose=verbose,
            criteria="normal",
            output_filename=os.path.join(
                compliance_dir,
                f"{os.path.basename(out_ncfile)}.compliance.html",
            ),
            output_format="html",
        )

    if plot:
        plot_dir = os.path.join(output_dir, "plots")
        create_dir(plot_dir)
        qc_plots.generate_plots(
            ncfile.filename,
            plot_dir,
            verbose,
        )

    logger.info("=" * 66)
    logger.info(
        "Done at %s %s"
        % (
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"),
            datetime.datetime.now().astimezone().tzinfo,
        ),
    )


if __name__ == "__main__":
    parser = clparser()
    args = parser.parse_args()
    # Positional arguments
    config = args.config
    input_file = args.input_file
    output_dir = args.output_dir
    # Optional
    plot = args.plot
    log = args.log
    verbose = args.verbose
    compliance_check = args.compliance_check

    if verbose:
        logger.addHandler(logging.StreamHandler())
    if log:
        if os.path.dirname(log) and not os.path.exists(os.path.dirname(log)):
            try:
                os.makedirs(os.path.dirname(log))
            except OSError:
                raise
        fh = logging.FileHandler(log)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    asv_ctd_qa(config, input_file, output_dir, plot, verbose, log, compliance_check)
