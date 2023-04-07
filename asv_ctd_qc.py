#!/usr/bin/env python

import argparse
import calendar
import csv
import datetime
import json
import logging
import os
import sys
import typing

import netCDF4 as nc  # noqa N813
import numpy as np
import numpy.typing as npt
import pandas as pd

import utils.alert as alert
import utils.ncdump as ncdump
import utils.qc_plots as qc_plots
from ioos_qc import qartod
from ioos_qc.config import QcConfig
from ioos_qc.qartod import QartodFlags
from ioos_qc.stores import PandasStore
from ioos_qc.streams import PandasStream

# from ioos_qc.utils import check_timestamps
from utils.compliance import Compliance

# Keyword to look for indicating the start of the data column row
INIT_COLUMN_ROW = "Date / Time"
# Log section separator
LOGGER_SEPARATOR = "=" * 80


class UnicodeReader:
    """A CSV reader which will iterate over lines in the CSV
    file "f" which is encoded in the given encoding.

    Args:
    ----
        f (str): file reader
        dialect (str, optional): File dialect. Defaults to csv.excel.
        encoding (str, optional): File encoding. Defaults to "utf-8".
        kwds (dict, optional): Optional keyword arguments.
    """

    def __init__(
        self: "UnicodeReader",
        f: str,
        dialect=csv.excel,  # noqa ANN001
        encoding: str = "utf-8",
        **kwds: dict,
    ) -> None:
        self.reader = csv.reader(f, dialect=dialect, **kwds)
        self.rows_read = 0

    def next(self: "UnicodeReader") -> list:
        """Return next reader object.

        Returns:
        -------
            (list): Row object.
        """
        return next(self.reader)

    def __iter__(self: "UnicodeReader") -> typing.Iterator:
        """Create an iterator from the provided object.

        Returns:
            iterator (iter): Iterator object.
        """
        return iter(self)


class FileParser:
    """Parse input file attributes.

    Args:
    ----
        fname (str): Input file name.
        header_rows (int): Number of input file header rows.
        column_delimiter (str): Input file column delimiter.
    """

    def __init__(
        self: "FileParser",
        fname: str,
        header_rows: int,
        column_delimiter: str | None = None,
    ) -> None:
        self.filename = fname
        self.delimiter = column_delimiter
        self.header_rows = header_rows
        self.datarows_start = self.header_rows + 1
        self.cols = []
        self.filerows = 0
        self.datarows = 0
        self.encoding = None
        self.reader = self.file_reader()
        self.eval_file()

    def file_reader(self: "FileParser"):  # noqa ANN201
        """Determine file dialect and encoding, then create.

        and return a reader object
        """
        self.dialect = csv.get_dialect("unix")
        try:
            self.dialect = csv.Sniffer().sniff(open(self.filename).readline())
        except Exception:
            raise
        try:
            if self.encoding:
                reader = UnicodeReader(
                    open(self.filename),
                    self.dialect,
                    self.encoding,
                ).reader
            else:
                reader = csv.reader(open(self.filename), self.dialect)
        except Exception:
            raise
        if self.header_rows > 0:
            try:
                for _row in range(self.header_rows):
                    next(reader)
            except Exception:
                raise
        self.cols = next(reader)
        return reader

    def eval_file(self: "FileParser") -> "FileParser":
        """Determine number of total rows and data rows in file."""
        self.filerows = self.header_rows + 1
        for datarow in self.reader:
            self.filerows += 1
            if len(datarow) > 0:
                self.datarows += 1
        return self

    def filestats(self: "FileParser") -> str:
        """Returns multi-line string of file information."""
        fstats = (
            """Filename: %s\nFile rows: %d\nHeader rows: %d\nColumn row: %d\nData rows: %d\n"""  # noqa
            % (
                self.filename,
                self.filerows,
                self.header_rows,
                self.datarows_start,
                self.datarows,
            )
        )
        fstats += "Columns:\n"
        for idx, col in enumerate(self.cols):
            i = f"0{idx}" if idx < 10 else str(idx)
            if col == self.cols[-1]:
                fstats += f"  {i} : {col}"
            else:
                fstats += f"  {i} : {col}\n"
        return fstats


class Parameters(FileParser):
    """Parameter parser.

    Args:
    ----
        config (dict): Global configuration.
        fname (str): Parameter file name.
        header_rows (int): Number of header rows in the parameter file
        preceeding the data column row.
        column_delimiter (str): Column delimiter of the parameter file.
    """

    def __init__(
        self: "Parameters",
        config: dict,
        input_file: str,
        num_header_rows: int,
        column_delimiter: str | None = None,
    ) -> None:
        super().__init__(input_file, num_header_rows, column_delimiter)
        self.config = config
        self.reader = self.file_reader()
        self.data_array = []

    def parameter_dataframe(self: "Parameters") -> None:
        """Convert parameters to dataframe."""
        for _idx, row in enumerate(self.reader):
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
        log_variables(keys)
        self.df.rename(columns=keys, inplace=True)
        self.df["syntax_test"] = self.syntax_test()
        self.df["gap_test"] = self.gap_test(
            pd.to_datetime(self.df["time"]),
            self.config["ioos_qc"]["global"]["gap_test"]["max_time_interval"],
        )
        for key in keys:
            self.df[keys[key]] = self.df[keys[key]].astype(
                self.config["parameters"][key]["dtype"],
            )

    def syntax_test(self: "Parameters") -> npt.ArrayLike:
        """Syntax test.

        Checks that data row is read as expected, matching the number
        of column headers in the data file.

        Returns
        -------
            np.typing.ArrayLike: Numpy array.
        """
        # Start with everything as passing (1)
        flag_arr = np.ma.ones(self.datarows, dtype="uint8")
        for idx, row in enumerate(self.file_reader()):
            row = tuple(filter(None, row))
            if len(row) != len(self.cols):
                flag_arr[idx] = QartodFlags.FAIL
        if np.all(flag_arr == 4):
            raise Exception("All syntax tests have failed.")
        return flag_arr

    def gap_test(
        self: "Parameters",
        times: np.ndarray,
        max_time_interval: int = 1,
    ) -> npt.ArrayLike:
        """Sanity checks for timestamp arrays.
        Checks that the times supplied are in monotonically increasing
        chronological order, and optionally that time intervals between
        measurements do not exceed a value `max_time_interval`.  Note that this is
        not a QARTOD test, but rather a utility test to make sure times are in the
        proper order and optionally do not have large gaps prior to processing the
        data.
        Args:
            times: Input array of timestamps
            max_time_interval: The interval between values should not exceed this
                value (seconds). Defaults to 1
        """
        flag_arr = np.ma.ones(len(times), dtype="uint8")
        fail_interval = np.timedelta64(max_time_interval, "s")

        time_diff = np.diff(times)
        if np.any(time_diff > fail_interval):
            for idx, time in enumerate(time_diff):
                if time > fail_interval:
                    flag_arr[idx + 1] = QartodFlags.FAIL
        return flag_arr


class NetCDF:
    """Create a NetCDF object and provide methods for adding
    self-describing attributes.
    """

    def __init__(
        self: "NetCDF",
        config: dict,
        outfile: str,
        data_model: str = "NETCDF4",
    ) -> None:
        """
        Args:
        ----
            config (dict): QC configuration json/dict.
            outfile (str): Output NetCDF file.
            data_model (str, optional): Data model. Defaults to "NETCDF4".
        """
        self.filename = outfile
        self.data_model = data_model
        self.rootgrp = nc.Dataset(self.filename, "w", format=self.data_model)
        self.config = config
        self.metadata()
        self.create_dimensions()

    def metadata(self: "NetCDF") -> None:
        """_summary_.

        Args:
        ----
            self (NetCDF): _description_
        """
        for attr in self.config["attributes"]:
            if attr == "history":
                setattr(
                    self.rootgrp,
                    attr,
                    f"""{self.config['attributes'][attr]} {
                        datetime.datetime.now().isoformat()}""",
                )
            else:
                setattr(self.rootgrp, attr, self.config["attributes"][attr])

    def add_metadata(self: "NetCDF", **kwargs: dict) -> None:
        if len(kwargs) > 0:
            for arg in kwargs:
                setattr(self.rootgrp, arg, kwargs[arg])

    def create_dimensions(self: "NetCDF") -> None:
        """Initialize the dataset dimensions."""
        self.rootgrp.createDimension("time", None)

    def create_ancillary_variables(
        self: "NetCDF",
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
            v[:] = data.to_numpy().astype("datetime64[s]")
        else:
            v[:] = data.to_numpy()

    def create_test_variables(
        self: "NetCDF",
        data: pd.DataFrame,
        ancillary_variable: str,
        test_name: str,
        variables: dict,
    ) -> None:
        """Write the varaible test results.

        Args:
        ----
            data (pd.DataFrame): Observations.
            ancillary_variable (str): Base parameters.
            test_name (str): Test name.
            variables (dict): Test variables.
        """
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
            elif type(variables[attr]) == pd.Series:  # noqa RET507
                setattr(v, attr, variables[attr.to_numpy()])
            elif attr == "config":
                for member in variables[attr].members:
                    v.config_fspan_minv = float(member.fspan.minv)
                    v.config_fspan_maxv = float(member.fspan.maxv)
                    v.config_tspan_minv = float(
                        (member.tspan.minv - pd.Timestamp("1970-01-01"))
                        // pd.Timedelta("1s"),
                    )
                    v.config_tspan_maxv = float(
                        (member.tspan.maxv - pd.Timestamp("1970-01-01"))
                        // pd.Timedelta("1s"),
                    )
                    v.config_vspan_minv = float(member.vspan.minv)
                    v.config_vspan_maxv = float(member.vspan.maxv)
                    v.config_zspan_minv = float(member.zspan.minv)
                    v.config_zspan_maxv = float(member.zspan.maxv)
            else:
                setattr(v, attr, variables[attr])

        if "syntax_test" in test_name:
            v[:] = data["syntax_test"].to_numpy()
        elif "gap_test" in test_name:
            v[:] = data["gap_test"].to_numpy()
        else:
            v[:] = data[test_name].to_numpy()


def create_dir(dir: str) -> None:
    """Execute the makedirs method on a given path.

    Arguments:
    ---------
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
        help="""Create an HTML file containing plots of QC flags.
        Files are stored under a subdirectory of the specified output_dir.""",
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


def log_variables(dict: dict) -> None:
    """Write keyword variables to the logger.

    Each key:value pair is written to a new line.
    """
    # Get longest key name and pad every other key with spaces for pretty formatting.
    max_chars = max(map(len, dict)) + 1
    if len(dict) > 0:
        for key in dict:
            n = max_chars - len(key)
            padding = " " * n
            logger.info(f"{key}{padding}: {dict[key]}")


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
                f"{today.year} 2 {eom} 23:59:59",
                "%Y %m %d %H:%M:%S",
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
    """Run QC tests.

    Args:
        config (dict): Configuration.
        parameters (Parameters): Parameters.
        ncfile (NetCDF): NetCDF file.
    """

    ioos_qc_config = config["ioos_qc"]

    for variable in ioos_qc_config["variables"]:
        variable_config = ioos_qc_config["variables"][variable]
        if "valid_range_test" in variable_config["axds"]:
            variable_config["axds"]["valid_range_test"]["valid_span"] = tuple(
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
            ncfile.add_metadata(
                season=season,
                window_start=window_start.isoformat(),
                window_end=window_end.isoformat(),
            )

            cc = qartod.ClimatologyConfig()
            cc.add(
                tspan=(window_start, window_end),
                fspan=tuple(
                    np.array(
                        variable_config["qartod"]["gross_range_test"]["fail_span"],
                    ),
                ),
                vspan=tuple(
                    np.array(
                        variable_config["qartod"]["gross_range_test"]["suspect_span"],
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

    # Run the tests and store results
    ps = PandasStream(parameters.df)
    qc = QcConfig(ioos_qc_config["variables"])
    qc_results = ps.run(qc)
    store = PandasStore(qc_results)
    results_store = store.save(write_data=False, write_axes=False)

    # Write results to netcdf
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

        # Add syntax test results. This isnt run by ioos_qc.
        var = f"{variable}_qartod_syntax_test"
        ncfile.create_test_variables(
            parameters.df,
            variable,
            f"{variable}_qartod_syntax_test",
            {},
        )
        results_store[var] = parameters.df["syntax_test"]

        # Add timestamp test results. This isnt run by ioos_qc.
        var = f"{variable}_qartod_gap_test"
        ncfile.create_test_variables(
            parameters.df,
            variable,
            f"{variable}_qartod_gap_test",
            config["ioos_qc"]["global"]["gap_test"],
        )
        results_store[var] = parameters.df["gap_test"]

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


def get_num_header_rows(input_file: str) -> int:
    """Open the input file and check for header rows.

    Searches the first item in each line for a substring that
    exactly matches  INIT_COLUMN_ROW. If the end of the file
    is reached without finding a match, a StopIteration error
    is raised and indicates how many rows were searched while not
    finding an exact string match.

    Args
    ----
        input_file (str): Input file to checck for header rows.

    Returns
    -------
        int: Number of header rows.
    """
    try:
        dialect = csv.Sniffer().sniff(open(input_file).readline())
    except Exception:
        raise
    try:
        reader = csv.reader(open(input_file), dialect)
    except Exception:
        raise
    num_header_rows = 0
    try:
        while True:
            row = next(reader)
            if row[0] == INIT_COLUMN_ROW:
                break
            num_header_rows += 1
    except StopIteration as error:
        raise Exception(
            """Could not determine start of data stream.
            Searched %i lines for string starting with %s."""
            % (num_header_rows, INIT_COLUMN_ROW),
        ) from error
    return num_header_rows


def asv_ctd_qc(
    config_file: str,
    input_file: str,
    output_dir: str,
    plot: bool = False,
    verbose: bool = False,
    compliance_check: bool = False,
) -> None:
    """Run IOOS QC tests.

    Args:
    ----
        - config (str): Configuration JSON file path.
        - input_file (str): Input file containing CTD parameters.
        - output_dir (str): Directory to save output files.
        - plot (bool, optional): Produce plots of CTD parameters with
          QARTOD flags in HTML format. Defaults to False.
        - verbose (bool, optional): Verbose mode. Defaults to False.
        - log (str, optional): Log file path. Defaults to None.
        - compliance_check (bool, optional): Run IOOS compliance checks
          on the output NetCDF file. Defaults to False.
    """
    # Validate input arguments
    if not os.path.exists(config_file):
        raise OSError("The input configuration file does not exist: %s", config_file)
    if not os.path.splitext(config_file)[1] == ".json":
        raise OSError("The input configuration file is invalid: %s" % config_file)
    if not os.path.exists(input_file):
        raise OSError("The input data file does not exist: %s" % input_file)

    # Initialize some variables
    input_file_basefn = os.path.basename(input_file)
    num_header_rows = get_num_header_rows(input_file)
    out_ncfile = os.path.join(output_dir, f"{input_file_basefn}.nc")

    # Write a log preamble
    logger.info(LOGGER_SEPARATOR)
    log_variables(
        {
            "script_start_time": "%s %s"
            % (
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"),
                datetime.datetime.now().astimezone().tzinfo,
            ),
            "config_file": config_file,
            "input_file": input_file,
            "input_file_basefn": input_file_basefn,
            "num_headerrows": num_header_rows,
            "output_dir": output_dir,
            "out_ncfile": out_ncfile,
            "plot": plot,
            "verbose": verbose,
        },
    )

    # Load JSON configuration
    try:
        with open(config_file, encoding="utf-8") as f:
            config_json = json.load(f)
    except OSError:
        raise

    logger.info(LOGGER_SEPARATOR)
    # Evaulate input parameters and create Parameter object.
    logger.info("EVAULATING INPUT FILE")
    parameters = Parameters(
        config=config_json,
        input_file=input_file,
        num_header_rows=num_header_rows,
    )
    # Log parameter file details
    logger.info(parameters.filestats())

    logger.info(LOGGER_SEPARATOR)
    # Create dataframe of input parameters
    # Dataframe will be used to append QARTOD flag columns
    parameters.parameter_dataframe()

    # Create directory to store output, if it doesnt exist already
    create_dir(output_dir)

    logger.info(LOGGER_SEPARATOR)
    # Initialize NetCDF file and write parameter data
    ncfile = NetCDF(config_json, out_ncfile)
    for parameter in config_json["parameters"]:
        logger.info(
            f"""CREATING ATTRIBUTES FOR {config_json['parameters'][
                parameter]['standard_name']}""",
        )
        ncfile.create_ancillary_variables(
            parameters.df[config_json["parameters"][parameter]["standard_name"]],
            parameter,
            config_json["parameters"][parameter],
        )

    logger.info(LOGGER_SEPARATOR)
    # Run IOOS_QC tests
    run_qc(config_json, parameters, ncfile)

    # Pretty print NetCDF file information, also store sections as variables
    nc_attrs, nc_dims, nc_vars = ncdump.ncdump(ncfile)

    if compliance_check:
        logger.info(LOGGER_SEPARATOR)
        logger.info("RUNNING CF COMPLIANCE CHECKS")

        compliance_dir = os.path.join(output_dir, "compliance_checks")
        output_filename = os.path.join(
            compliance_dir,
            f"{os.path.basename(out_ncfile)}.compliance.html",
        )
        create_dir(compliance_dir)
        convention = json.dumps(config_json["attributes"]["Conventions"]).strip('"')
        Compliance(convention=convention).run_checker(
            out_ncfile,
            output_filename,
            verbose,
        )
        logger.info(f"Compliance results saved to: {output_filename}")

    if plot:
        logger.info(LOGGER_SEPARATOR)
        plot_dir = os.path.join(output_dir, "plots")
        create_dir(plot_dir)
        qc_plots.generate_plots(
            ncfile.filename,
            plot_dir,
        )
        logger.info(f"Plots saved to : {plot_dir}")

    logger.info(LOGGER_SEPARATOR)
    log_variables(
        {
            "Done at": "%s %s"
            % (
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"),
                datetime.datetime.now().astimezone().tzinfo,
            ),
        },
    )

    # Alert staff on successful data exchange or critical issues.
    alert.main(input_file=input_file)


if __name__ == "__main__":
    parser = clparser()
    args = parser.parse_args()
    # Positional arguments
    config_file = args.config
    input_file = args.input_file
    output_dir = args.output_dir
    # Optional
    plot = args.plot
    log = args.log
    verbose = args.verbose
    compliance_check = args.compliance_check

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s : %(msecs)04d : %(name)s : %(levelname)s : %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

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

    try:
        asv_ctd_qc(
            config_file=config_file,
            input_file=input_file,
            output_dir=output_dir,
            plot=plot,
            verbose=verbose,
            compliance_check=compliance_check,
        )
    except Exception:
        alert.send_alert(
            subject="ASV CTD QC Alert",
            body=f"QC workflow failed to run on input file {input_file}.",
        )
        sys.exit(1)
