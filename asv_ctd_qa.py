#!/usr/bin/env python
# coding=utf-8
import argparse
import calendar
import csv
import datetime
import json
import logging
import os
import sys

import netCDF4 as nc
import numpy as np
import pandas as pd
from cf_units import Unit  # noqa

import ncdump
import qc_plots
from ioos_qc import qartod
from ioos_qc.config import QcConfig
from ioos_qc.qartod import QartodFlags
from ioos_qc.stores import PandasStore
from ioos_qc.streams import PandasStream

ERROR_LEVEL = {0: "INFO", 1: "WARNING", 2: "CRITICAL"}
error_list = []

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(asctime)s : %(msecs)04d : %(name)s : %(levelname)s : %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log_dir = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "logs")
if not os.path.exists(log_dir):
    try:
        os.makedirs(log_dir)
    except OSError:
        raise
log_file = f"log_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
log_path = os.path.join(log_dir, log_file)
fh = logging.FileHandler(log_path)
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
logger.addHandler(fh)


def print_errors(log: logging.Logger, error_list: list) -> None:
    error_str = ""
    for idx, error in enumerate(error_list):
        error_str = error_str + f"{idx}.  {error[0]}: {error[1]}\n"
    log.warning(error_str)


class ScriptError(Exception):
    """Script error tracker."""

    def __init__(self, errlvl, errmsg, filename=None):
        self.errlvl = errlvl
        self.errmsg = errmsg
        self.filename = filename
        error_list.append((self.errlvl, self.errmsg, self.filename))


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
            ScriptError(
                ERROR_LEVEL[2],
                "Can't open file to determine format.",
                self.filename,
            )
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
            ScriptError(ERROR_LEVEL[2], "Can't open file to read data.", self.filename)
        if self.headrows > 0:
            try:
                for row in range(self.headrows):
                    next(reader)
            except Exception:
                ScriptError(
                    ERROR_LEVEL[2],
                    "Can't read column header line.",
                    self.filename,
                )
        self.cols = next(reader)
        # try:
        #     # Remove blank columns that are result of trailing delimiters
        #     self.cols.remove("")
        # except Exception:
        #     pass
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
            logger.error("Unable to parse input data.")
            raise ScriptError(ERROR_LEVEL[2], "Unable to parse input data.")

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
            raise ScriptError(ERROR_LEVEL[2], "All syntax tests have failed.")
        return flag_arr


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
            dtype = "int64"
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

        # Include syntax_test, rollup_qc, location_test, check_timestamps

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
                # for c in attr:
                # print(c, attr[c])
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
            raise ScriptError(ERROR_LEVEL[2], "Cannot create dir.", dir)


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
        "header_rows",
        help="""Number of rows preceeding the row
                        containing column headers.""",
    )
    parser.add_argument(
        "output_dir",
        help="Directory for output files.",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        dest="plot",
        default=False,
        help="Create an HTML file containing plots of QC flags.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="Control the amount of information to display.",
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

            # Remove this
            # ================
            # window_start = datetime.datetime.strptime("2020-01-01", "%Y-%m-%d")
            # window_end = datetime.datetime.strptime("2022-12-31", "%Y-%m-%d")
            # ================

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
            # print(cc.members)
            variable_config["qartod"]["climatology_test"]["config"] = cc
        # Add global tests
        variable_config["qartod"].update(config["ioos_qc"]["global"]["qartod"])

    ps = PandasStream(parameters.df)
    qc = QcConfig(ioos_qc_config["variables"])
    qc_results = ps.run(qc)
    store = PandasStore(qc_results)
    store.compute_aggregate()
    results_store = store.save(write_data=False, write_axes=False)

    for variable in ioos_qc_config["variables"]:
        for test_type in ioos_qc_config["variables"][variable]:
            for test in ioos_qc_config["variables"][variable][test_type]:
                test_name = f"{variable}_{test_type}_{test}"
                # print(ioos_qc_config["variables"][variable][test_type])
                logger.info(f"RUNNING TEST {test_name}")
                if ioos_qc_config["variables"][variable][test_type][test] is not None:
                    ncfile.create_test_variables(
                        results_store,
                        variable,
                        test_name,
                        ioos_qc_config["variables"][variable][test_type][test],
                    )


def asv_ctd_qa(
    config: str,
    input_file: str,
    num_headerrows: int,
    output_dir: str,
    plot: bool = False,
    verbose: bool = False,
) -> None:
    """Run IOOS QC tests."""

    # Validate input arguments
    if not os.path.exists(config):
        raise ScriptError(
            ERROR_LEVEL[2],
            "The input configuration file does not exist.",
            config,
        )
    if not os.path.splitext(config)[1] == ".json":
        raise ScriptError(
            ERROR_LEVEL[2],
            "The input configuration file is invalid",
            config,
        )
    if not os.path.exists(input_file):
        raise ScriptError(
            ERROR_LEVEL[2],
            "The input data file does not exist.",
            input_file,
        )
    try:
        num_headerrows = int(num_headerrows)
    except ValueError:
        raise ScriptError(ERROR_LEVEL[2], "The [num_headerrows] argument is invalid.")

    input_file_basefn = os.path.basename(input_file)
    create_dir(output_dir)
    ncfile_dir = os.path.join(output_dir, "netcdf")
    create_dir(ncfile_dir)
    out_ncfile = os.path.join(ncfile_dir, f"{input_file_basefn}.nc")

    if verbose:
        logger.addHandler(logging.StreamHandler())
    preamble(
        config=config,
        input_file=input_file,
        input_file_basefn=input_file_basefn,
        num_headerrows=num_headerrows,
        output_dir=output_dir,
        ncfile_dir=ncfile_dir,
        out_ncfile=out_ncfile,
        log_dir=log_dir,
        plot=plot,
        verbose=verbose,
    )

    try:
        with open(config, "r", encoding="utf-8") as f:
            config_json = json.load(f)
    except OSError:
        raise ScriptError(ERROR_LEVEL[2], "Could not read config file.", config)

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

    ncfile = NetCDF(config_json, out_ncfile)
    for parameter in config_json["parameters"]:
        logger.info(f"CREATING ATTRIBUTES FOR {parameter}")
        ncfile.create_ancillary_variables(
            parameters.df[config_json["parameters"][parameter]["standard_name"]],
            parameter,
            config_json["parameters"][parameter],
        )

    run_qc(config_json, parameters, ncfile)

    # Pretty print NetCDF file information, also store sections as variables
    logger.info("=" * 66)
    logger.info("LOGGING OUTPUT NETCDF ATTRIBUTES")
    nc_attrs, nc_dims, nc_vars = ncdump.ncdump(ncfile, verb=True if verbose else False)

    if plot:
        logger.info("=" * 66)
        logger.info("CREATING OBSERVATION PLOTS WITH QC FLAGS")
        plot_dir = os.path.join(output_dir, "plots")
        create_dir(plot_dir)
        qc_plots.generate_plots(
            ncfile.filename,
            plot_dir,
            verbose,
        )

    # Check for errors to report on
    if len(error_list) > 0:
        if verbose:
            print_errors(error_list)

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
    num_headerrows = args.header_rows
    output_dir = args.output_dir
    # Optional
    plot = args.plot
    verbose = args.verbose
    asv_ctd_qa(config, input_file, num_headerrows, output_dir, plot, verbose)
