#! /usr/bin/python
# -*- coding: utf-8 -*-
import argparse
import csv
import datetime
import logging
import os
import smtplib
import socket
import sqlite3

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import quantities as pq
import xarray as xr
import xlsx2csv as xl

from ioos_qc import argo, axds
from ioos_qc import qartod as qc
from ioos_qc import utils

NOTIFIER_EMAIL_FROM = "cgrant@integral-corp.com"
NOTIFIER_EMAIL_TO = ["cgrant@integral-corp.com"]
NOTIFIER_EMAIL_CC = []
HOST = socket.gethostbyname(socket.gethostname())

ERROR_LEVEL = {0: "INFO", 1: "WARNING", 2: "CRITICAL"}
error_list = []


def print_errors(log, error_list):
    error_str = ""
    for idx, error in enumerate(error_list):
        error_str = error_str + f"{idx}.  {error[0]}: {error[1]}\n"
    log.warning(error_str)


class ScriptError(Exception):
    """Script error tracker"""

    def __init__(self, errlvl, errmsg, filename=None):
        self.errlvl = errlvl
        self.errmsg = errmsg
        self.filename = filename
        error_list.append((self.errlvl, self.errmsg, self.filename))


class Notifier:
    """Send SMTP notification"""

    def __init__(
        self, _host: str, _from: str, _to: list, _cc: list, subj: str, msg: str
    ) -> None:
        self._from = _from
        self._to = _to
        self._cc = _cc
        self.subj = subj
        self._host = _host
        self.message = f"""
            From: {self._from}\n
            To: {self._to}\n
            CC: {self._cc}\n
            Subject: {self.subj}\n\n
            {msg}"""

    def send(self):
        try:
            with smtplib.SMTP(self._host, 1025) as server:
                server.sendmail(self._from, self._to + self._cc, self.message)
        except Exception:
            raise ScriptError(
                ERROR_LEVEL[1], "Failed to send notification email using SMTP."
            )


def clparser():
    """Create a parser to handle input arguments and displaying
    a script specific help message."""
    desc_msg = """Evaulate a data file of ASV CTD readings
        and applies quality assurence
        checks following QARTOD methods and assigning data
        quality flags as appropriate. Transform
        results into NetCDF format following IC standards."""
    parser = argparse.ArgumentParser(description=desc_msg)
    parser.add_argument("input_file", help="Path to the input sensor data file.")
    parser.add_argument(
        "header_rows",
        help="""Number of rows preceeding the row
                        containing column headers.""",
    )
    parser.add_argument("param_file", help="Path to sensor threshold parameters file.")
    parser.add_argument("output_dir", help="Path for output files to be stored.")
    parser.add_argument("log_dir", help="Directory to store log files.")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="Control the amount of information to display.",
    )
    parser.add_argument(
        "-e",
        "--error_monitoring",
        action="store_true",
        dest="error_monitoring",
        default=False,
        help="""Treat errors as CRITICAL and notify operators of
                any issues via email. This should generally only
                be used in production.""",
    )
    return parser


def ncdump(nc_fid, log, verb=True):
    """
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    """

    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            log.info("\t\ttype: %s" % repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                log.info(
                    "\t\t%s: %s"
                    % (ncattr, repr(nc_fid.variables[key].getncattr(ncattr)))
                )
        except KeyError:
            log.info("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        log.info("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            log.info("\t%s: %s" % (nc_attr, repr(nc_fid.getncattr(nc_attr))))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        log.info("NetCDF dimension information:")
        for dim in nc_dims:
            log.info("\tName: %s" % dim)
            log.info("\t\tsize: %s" % len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        log.info("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                log.info("\tName: %s" % var)
                log.info("\t\tdimensions: %s" % nc_fid.variables[var].dimensions)
                log.info("\t\tsize: %s" % nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars


def create_dirs(path: str) -> None:
    """Execute the makedirs method on a given path"""
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except Exception as err:
            raise err


def start_logging(log_dir) -> logging.Logger:
    """Create a logger instance.

    Arguments:
      dir: Directory to store log files.
    """
    create_dirs(log_dir)
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s : %(msecs)04d : %(levelname)s : %(message)s",
        filename=(
            os.path.join(
                log_dir, f"log_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
            )
        ),
        filemode="w",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger(__name__)


def log_preamble(log: logging.Logger, **kwargs) -> None:
    """Write a preamble to a log file consisting of key:value pairs."""
    log.info("=" * 75)
    log.info(
        "start time: %s %s"
        % (
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"),
            datetime.datetime.now().astimezone().tzinfo,
        )
    )
    if len(kwargs) > 0:
        for arg in kwargs:
            log.info("%s: %s" % (arg, kwargs[arg]))


class UnicodeReader:
    """A CSV reader which will iterate over lines in the CSV file "f",
    which is encoded in the given encoding."""

    def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
        self.reader = csv.reader(f, dialect=dialect, **kwds)
        self.rows_read = 0

    def next(self):
        row = next(self.reader)
        return row

    def __iter__(self):
        return iter(self)


class FileParser:
    """Parse input file attributes"""

    def __init__(self, fname: str, header_rows: int, column_delimiter: str) -> None:
        self.filename = fname
        self.delimiter = column_delimiter
        self.headrows = header_rows
        self.colrow = self.headrows + 1
        self.cols = []
        self.filerows = 0
        self.datarows = 0
        self.encoding = None
        self.reader = self.file_reader()
        self.eval_file()

    def file_reader(self):
        """Determine file dialect and encoding, then create
        and return a reader object"""
        try:
            dialect = csv.Sniffer().sniff(open(self.filename, "rt").readline())
        except Exception:
            ScriptError(
                ERROR_LEVEL[2], "Can't open file to determine format.", self.filename
            )
        try:
            if self.encoding:
                reader = UnicodeReader(
                    open(self.filename, "rt"), dialect, self.encoding
                )
            else:
                reader = csv.reader(open(self.filename, "rt"), dialect)
        except Exception:
            ScriptError(ERROR_LEVEL[2], "Can't open file to read data.", self.filename)
        if self.headrows > 0:
            try:
                for row in range(self.headrows):
                    next(reader)
            except Exception:
                ScriptError(
                    ERROR_LEVEL[2], "Can't read column header line.", self.filename
                )
        self.cols = next(reader)
        try:
            # Remove blank columns that are result of trailing delimiters
            self.cols.remove("")
        except Exception:
            pass
        return reader

    def eval_file(self):
        """Determine number of total rows and data rows in file"""
        self.filerows = self.headrows + 1
        for datarow in self.reader:
            self.filerows += 1
            if len(datarow) > 0:
                self.datarows += 1
        return self

    def filestats(self):
        """Returns multi-line string of file information"""
        fstats = (
            "Filename: %s\nFile rows: %d\nHeader rows: %d\nColumn row: %d\nData rows: %d\n"
            % (self.filename, self.filerows, self.headrows, self.colrow, self.datarows)
        )
        fstats += "Columns:\n"
        for col in self.cols:
            fstats += "  %s\n" % col
        return fstats


class MakeDataFrame:
    def __init__(self, fileobj):
        self.fobj = fileobj
        self.cols = fileobj.cols
        self.reader = fileobj.file_reader()

    def sensor_dataframe(self):
        """Create a Pandas DataFrame from an array of arrays retrieved from the reader"""
        self.data_array = []
        self.syntax_test_failed_rows = []

        for (idx, row) in enumerate(self.reader):
            row = list(filter(None, row))
            if len(row) < len(self.cols):
                if idx + self.fobj.colrow + 1 not in self.syntax_test_failed_rows:
                    self.syntax_test_failed_rows.append(idx + self.fobj.colrow + 1)
                ScriptError(
                    ERROR_LEVEL[1],
                    f"""The input file is missing {len(self.cols) - len(row)} columns on line {idx + self.fobj.colrow + 1}.""",
                    self.fobj.filename,
                )
                for c in range(len(self.cols) - len(row)):
                    row.append(np.nan)
            if len(row) > len(self.cols):
                if idx + self.fobj.colrow + 1 not in self.syntax_test_failed_rows:
                    self.syntax_test_failed_rows.append(idx + self.fobj.colrow + 1)
                ScriptError(
                    ERROR_LEVEL[1],
                    f"""The input file is missing {len(row) - len(self.cols)} column headers on line {idx + self.fobj.colrow + 1}. Columns with no headers will be removed.""",
                    self.fobj.filename,
                )
                row = row[: -(len(row) - len(self.cols))]
            self.data_array.append(row)
        self.df = pd.DataFrame(self.data_array, columns=self.cols)
        """Set data types of the DataFrame"""
        for col in range(len(self.cols)):
            if col == 0:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype(
                    "datetime64"
                )
                self.df[self.df.columns[col]] = self.df[
                    self.df.columns[col]
                ] - datetime.datetime(2010, 1, 1, 0, 0)
                self.df[self.df.columns[col]] = self.df[
                    self.df.columns[col]
                ].dt.total_seconds()
            else:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype(
                    "float"
                )

    def parameter_dataframe(self):
        """Create a Pandas DataFrame from an array
        of arrays retrieved from the reader"""
        self.data_array = []
        for (idx, row) in enumerate(self.reader):
            self.data_array.append(list(filter(None, row)))
        self.df = pd.DataFrame(self.data_array, columns=self.cols)
        """Set data types of the DataFrame"""
        for col in range(len(self.cols)):
            if "value" in self.cols[col]:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype(
                    "float"
                )
            else:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype(
                    "str"
                )


def convert_xlsx(file: str) -> str:
    """Convert an XLSX file to CSV."""
    if not os.path.exists(file):
        raise ScriptError(ERROR_LEVEL[2], "The input file does not exist.", file)
    filename_noext = os.path.splitext(os.path.basename(os.path.abspath(file)))[0]
    output_csv = os.path.join(
        os.path.dirname(os.path.abspath(file)), f"{filename_noext}.csv"
    )
    xl.Xlsx2csv(file).convert(output_csv, sheetid=2)
    return output_csv


class NetCDF:
    """Create a NetCDF object and provide methods for adding self-describing attributes"""

    def __init__(self, outfile, data_model="NETCDF4"):
        self.filename = outfile
        self.data_model = data_model
        self.rootgrp = nc.Dataset(self.filename, "w", format=self.data_model)

    def metadata(self):
        self.rootgrp.title = "Hypoxia Monitoring Data"
        self.rootgrp.description = "NOAA ASV sensor parameter QA"
        self.rootgrp.history = "Created %s" % (
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        )
        self.rootgrp.source = "L3"
        self.rootgrp.conventions = "CF-1.6"
        self.rootgrp.author = "Integral Consulting Inc."
        # self.rootgrp.sensor_id = sensor_id
        # self.rootgrp.sensor_sn = sensor_sn

    def operational_variables(self, data):
        tim_dim = self.rootgrp.createDimension("time", None)
        tim = self.rootgrp.createVariable("time", np.dtype("float64").char, ("time",))
        tim.long_name = "Time (PST) as Seconds Since 2010-01-01 00:00:00"
        tim.standard_name = "time"
        tim.units = "seconds since 2010-01-01 00:00:00"
        tim.time_zone = "PST"
        tim.calendar = "standard"
        tim[:] = data.df["Date / Time"]

        lat = self.rootgrp.createVariable("lat", np.dtype("float64").char, ("time",))
        lat.long_name = "Latitude degrees"
        lat.standard_name = "latitude"
        lat.units = "degrees_north"
        lat[:] = data.df["LATITUDE;DEG"]

        lon = self.rootgrp.createVariable("lon", np.dtype("float64").char, ("time",))
        lon.long_name = "Longitude degrees"
        lon.standard_name = "longitude"
        lon.units = "degrees_east"
        lon[:] = data.df["LONGITUDE;DEG"]


def geophysical_variables():
    """
    Returns a dict of variable names matching the geophysical variable's
    standard names
    """
    standard_name_lookup = {
        "qartod": {
            "TEMPERATURE;C": "sea_water_temperature",
            "CONDUCTIVITY;MS/CM": "sea_water_electrical_conductivity",
            "PRESSURE;DBAR": "sea_water_pressure",
            "Calc. SALINITY; PSU": "sea_water_practical_salinity",
            "PH;PH": "sea_water_ph",
            "DISSOLVED OXYGEN;SAT%": "sea_water_dissolved_oxygen",
            "FLUOROMETER (C);UG/L": "sea_water_fluorescence",
            "TURBIDITY;FTU": "sea_water_turbidity",
            "ALTITUDE;M": "altitude",
            "Calc. DEPTH;M": "depth",
        },
        "operational": {
            "Date / Time": "datetime",
            "LATITUDE;DEG": "latitude",
            "LONGITUDE;DEG": "longitude",
        },
    }
    return standard_name_lookup


def syntax_test(times, failed_rows):
    flags = np.ma.ones(times.size, dtype="uint8")
    for i, v in enumerate(times):
        if i in failed_rows:
            flags[i] = qc.QartodFlags.FAIL
    return flags


class SensorQC(object):
    def __init__(self, data, config_params, ncfile, syntax_test_failed_rows):
        self.data = data
        self.params = config_params
        self.ncfile = ncfile
        self.syntax_test_failed_rows = syntax_test_failed_rows

    def create_variables(self, log, ncvariable, varname):
        """
        Returns a list of variable names for the newly created variables
        """
        name = ncvariable
        standard_name = ncvariable
        dims = "time"
        units = varname.split(";")[1].strip()
        source_name = varname
        log.info("Creating variables for %s", name)
        variables = []

        variable_name = name
        ncvar = self.ncfile.createVariable(
            variable_name, np.float64, dims, fill_value=np.int8(9)
        )
        ncvar.source_name = source_name
        ncvar.units = units
        ncvar.standard_name = standard_name
        ncvar.long_name = standard_name
        ncvar[:] = self.data.df[source_name]
        variables.append(variable_name)
        return variables

    def find_qc_flags(self, ncvariable):
        """
        Returns a list of non-GliderDAC QC flags associated with a variable

        :param netCDF4.Variable ncvariable: Variable to get the status flag
                                            variables for
        """
        valid_variables = []
        for varname in self.ncfile.variables:
            if hasattr(self.ncfile.variables[varname], "standard_name"):
                stdvar = self.ncfile.variables[varname].getncattr("standard_name")
                if isinstance(stdvar, str):
                    var = stdvar.split(" ")
                else:
                    var = []
                if (
                    ncvariable in var
                    and "status_flag" in var
                    and "primary_flag" not in varname
                ):
                    valid_variables.append(varname)
        return valid_variables

    def create_qc_variables(self, log, ncvariable, varname):
        """
        Returns a list of variable names for the newly created variables for QC flags
        """
        name = ncvariable
        standard_name = ncvariable
        dims = "time"
        log.info("Creating QARTOD variables for %s", name)

        templates = {
            "gap": {
                "name": "qartod_%(name)s_gap_flag",
                "long_name": "QARTOD Gap Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "gap",
                "comment": "ioos_qartod",
            },
            "valid_range": {
                "name": "qartod_%(name)s_valid_range_flag",
                "long_name": "QARTOD Valid Range Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "valid_range",
                "comment": "ioos_qartod",
            },
            "location": {
                "name": "qartod_%(name)s_location_flag",
                "long_name": "QARTOD Location Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "location",
                "comment": "ioos_qartod",
            },
            "syntax": {
                "name": "qartod_%(name)s_syntax_flag",
                "long_name": "QARTOD Syntax Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "syntax",
                "comment": "ioos_qartod",
            },
            "climatological": {
                "name": "qartod_%(name)s_climatological_flag",
                "long_name": "QARTOD Climatological Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "climatological",
                "comment": "ioos_qartod",
            },
            "flat_line": {
                "name": "qartod_%(name)s_flat_line_flag",
                "long_name": "QARTOD Flat Line Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "flat_line",
                "comment": "ioos_qartod",
            },
            "gross_range": {
                "name": "qartod_%(name)s_gross_range_flag",
                "long_name": "QARTOD Gross Range Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "gross_range",
                "comment": "ioos_qartod",
            },
            "rate_of_change": {
                "name": "qartod_%(name)s_rate_of_change_flag",
                "long_name": "QARTOD Rate of Change Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "rate_of_change",
                "comment": "ioos_qartod",
            },
            "spike": {
                "name": "qartod_%(name)s_spike_flag",
                "long_name": "QARTOD Spike Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "spike",
                "comment": "ioos_qartod",
            },
            "pressure": {
                "name": "qartod_%(name)s_pressure_flag",
                "long_name": "QARTOD Pressure Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "pressure",
                "comment": "ioos_qartod",
            },
            "attenuated_signal": {
                "name": "qartod_%(name)s_attenuated_signal_flag",
                "long_name": "QARTOD Attenuated Signal Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "attenuated_signal",
                "comment": "ioos_qartod",
            },
            "density_inversion": {
                "name": "qartod_%(name)s_density_inversion_flag",
                "long_name": "QARTOD Density Inversion Test for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "qartod_test": "density_inversion",
                "comment": "ioos_qartod",
            },
            "primary": {
                "name": "qartod_%(name)s_primary_flag",
                "long_name": "QARTOD Primary Flag for %(standard_name)s",
                "standard_name": "%(standard_name)s status_flag",
                "flag_values": np.array([1, 2, 3, 4, 9], dtype=np.int8),
                "flag_meanings": "GOOD NOT_EVALUATED SUSPECT BAD MISSING",
                "comment": "ioos_qartod_primary",
            },
        }

        qcvariables = []

        for tname, template in list(templates.items()):
            if tname == "pressure" and standard_name != "sea_water_pressure":
                continue
            variable_name = template["name"] % {"name": name}

            if variable_name not in geophysical_variables()["qartod"]:
                ncvar = self.ncfile.createVariable(
                    variable_name, np.int8, dims, fill_value=np.int8(9)
                )
            else:
                ncvar = self.ncfile.variables[variable_name]

            ncvar.units = "1"
            ncvar.standard_name = template["standard_name"] % {
                "standard_name": standard_name
            }
            ncvar.long_name = template["long_name"] % {"standard_name": standard_name}
            ncvar.flag_values = template["flag_values"]
            ncvar.flag_meanings = template["flag_meanings"]
            # ncvar.references = template['references']
            if "qartod_test" in template:
                ncvar.qartod_test = template["qartod_test"]
            qcvariables.append(variable_name)

        return qcvariables

    def apply_qc(self, syntax_test_failed_rows, ncvariable):
        """
        Applies QC to a qartod variable
        """
        qc_tests = {
            "gap": utils.check_timestamps,
            "valid_range": axds.valid_range_test,
            "location": qc.location_test,
            "syntax": syntax_test,
            "climatological": qc.climatology_test,
            "flat_line": qc.flat_line_test,
            "gross_range": qc.gross_range_test,
            "rate_of_change": qc.rate_of_change_test,
            "spike": qc.spike_test,
            "attenuated_signal": qc.attenuated_signal_test,
            "pressure": argo.pressure_increasing_test,
            "density_inversion": qc.density_inversion_test,
        }

        qartod_test = getattr(ncvariable, "qartod_test", None)
        if not qartod_test:
            return
        standard_name = getattr(ncvariable, "standard_name").split(" ")[0]
        parent = self.ncfile.get_variables_by_attributes(standard_name=standard_name)[0]
        times, values, mask = self.get_unmasked(parent)
        # There's no data to QC
        if len(values) == 0:
            return

        qa_config = self.params.df[
            self.params.df["sensor"] == getattr(parent, "source_name")
        ]
        op_config = self.params.df[self.params.df["sensor"] == "OPERATOR"]

        test_params = {}

        if qartod_test == "syntax":
            # Use syntax_test_failed_rows (list of row indexes) to
            # assign fail or pass flags to each row
            test_params["times"] = times
            test_params["failed_rows"] = syntax_test_failed_rows

        if qartod_test == "gap":
            test_params["times"] = ma.getdata(times[~mask])

        if qartod_test == "location":
            test_params["lat"] = self.ncfile.get_variables_by_attributes(
                standard_name="latitude"
            )[0][:]
            test_params["lon"] = self.ncfile.get_variables_by_attributes(
                standard_name="longitude"
            )[0][:]
            test_params["bbox"] = tuple(
                [
                    op_config.loc[
                        op_config["parameter"] == "long_min", "parameter_value"
                    ].iloc[0],
                    op_config.loc[
                        op_config["parameter"] == "lat_min", "parameter_value"
                    ].iloc[0],
                    op_config.loc[
                        op_config["parameter"] == "long_max", "parameter_value"
                    ].iloc[0],
                    op_config.loc[
                        op_config["parameter"] == "lat_max", "parameter_value"
                    ].iloc[0],
                ]
            )

        if qartod_test == "rate_of_change":
            # test_params["inp"] = values
            test_params["inp"] = ma.getdata(values[~mask])
            time_units = self.ncfile.variables["time"].units
            # dates = np.array(nc.num2date(times, time_units), dtype="datetime64[ms]")
            test_params["tinp"] = ma.getdata(times[~mask])
            n_dev = qa_config.loc[
                qa_config["parameter"] == "n_deviation", "parameter_value"
            ].iloc[0]
            test_params["threshold"] = self.get_rate_of_change_threshold(
                ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev
            )

        if qartod_test == "climatological":
            cc = qc.ClimatologyConfig()
            cc.add(
                # Seasonal window
                tspan=tuple(
                    [
                        np.array(
                            nc.num2date(times, self.ncfile.variables["time"].units),
                            dtype="datetime64[ms]",
                        ).min(),
                        np.array(
                            nc.num2date(times, self.ncfile.variables["time"].units),
                            dtype="datetime64[ms]",
                        ).max(),
                    ]
                ),
                # Valid value range
                vspan=tuple(
                    [
                        qa_config.loc[
                            qa_config["parameter"] == "climate_min", "parameter_value"
                        ].iloc[0],
                        qa_config.loc[
                            qa_config["parameter"] == "climate_max", "parameter_value"
                        ].iloc[0],
                    ]
                ),
            )
            test_params["config"] = cc
            test_params["tinp"] = ma.getdata(times[~mask])
            test_params["inp"] = ma.getdata(values[~mask])
            test_params["zinp"] = np.array(self.data.df["Calc. DEPTH;M"])

        if qartod_test == "spike":
            test_params["inp"] = ma.getdata(times[~mask])
            test_params["suspect_threshold"] = qa_config.loc[
                qa_config["parameter"] == "spike_low", "parameter_value"
            ].iloc[0]
            test_params["fail_threshold"] = qa_config.loc[
                qa_config["parameter"] == "spike_high", "parameter_value"
            ].iloc[0]

        if qartod_test == "gross_range":
            test_params["inp"] = ma.getdata(values[~mask])
            test_params["fail_span"] = tuple(
                [
                    qa_config.loc[
                        qa_config["parameter"] == "sensor_min", "parameter_value"
                    ].iloc[0],
                    qa_config.loc[
                        qa_config["parameter"] == "sensor_max", "parameter_value"
                    ].iloc[0],
                ]
            )
            test_params["suspect_span"] = tuple(
                [
                    qa_config.loc[
                        qa_config["parameter"] == "climate_min", "parameter_value"
                    ].iloc[0],
                    qa_config.loc[
                        qa_config["parameter"] == "climate_max", "parameter_value"
                    ].iloc[0],
                ]
            )

        if qartod_test == "valid_range":
            test_params["inp"] = ma.getdata(values[~mask])
            test_params["valid_span"] = tuple(
                [
                    qa_config.loc[
                        qa_config["parameter"] == "sensor_min", "parameter_value"
                    ].iloc[0],
                    qa_config.loc[
                        qa_config["parameter"] == "sensor_max", "parameter_value"
                    ].iloc[0],
                ]
            )

        if qartod_test == "flat_line":
            test_params["inp"] = ma.getdata(values[~mask])
            test_params["tinp"] = ma.getdata(times[~mask])
            test_params["suspect_threshold"] = qa_config.loc[
                qa_config["parameter"] == "rep_cnt_susp", "parameter_value"
            ].iloc[0]
            test_params["fail_threshold"] = qa_config.loc[
                qa_config["parameter"] == "rep_cnt_fail", "parameter_value"
            ].iloc[0]
            test_params["tolerance"] = qa_config.loc[
                qa_config["parameter"] == "eps", "parameter_value"
            ].iloc[0]

        if qartod_test == "pressure":
            test_params["inp"] = ma.getdata(values[~mask])

        if qartod_test == "density_inversion":
            test_params["inp"] = ma.getdata(values[~mask])
            test_params["zinp"] = self.data.df["PRESSURE;DBAR"].to_numpy()

        if qartod_test == "attenuated_signal":
            test_params["inp"] = ma.getdata(values[~mask])
            test_params["tinp"] = ma.getdata(times[~mask])
            test_params["suspect_threshold"] = qa_config.loc[
                qa_config["parameter"] == "climate_min", "parameter_value"
            ].iloc[0]
            test_params["fail_threshold"] = qa_config.loc[
                qa_config["parameter"] == "climate_max", "parameter_value"
            ].iloc[0]
            test_params["min_obs"] = qa_config.loc[
                qa_config["parameter"] == "rep_cnt_susp", "parameter_value"
            ].iloc[0]

        if qartod_test == "climatological":
            qc_flags = qc_tests[qartod_test](
                config=test_params["config"],
                inp=test_params["inp"],
                tinp=test_params["tinp"],
                zinp=test_params["zinp"],
            )
        else:
            qc_flags = qc_tests[qartod_test](**test_params)
        ncvariable[~mask] = qc_flags

        for test_param in test_params:
            # if test_param in ("inp", "tinp", "zinp"):
            # continue

            if qartod_test == "climatological" and test_param == "config":
                for member in test_params[test_param].members:
                    ncvariable.setncattr("config_tspan", str(member.tspan))
                    ncvariable.setncattr("config_fspan", str(member.fspan))
                    ncvariable.setncattr("config_vspan", str(member.vspan))
                    ncvariable.setncattr("config_zspan", str(member.zspan))
                    ncvariable.setncattr("config_period", str(member.period))
            else:
                ncvariable.setncattr(test_param, test_params[test_param])

    def get_rate_of_change_threshold(
        self, values, times, time_units="seconds since 1970-01-01T00:00:00Z", n_dev=3
    ):
        """
        Return the threshold used for the rate of change test

        :param values: numpy array of values
        :param times: numpy array of times
        :param time_units: string defining time units
        """
        n_dev = n_dev  # Set to 3 standard deviations
        std = np.nanstd(values)
        thresh = n_dev * std
        tstep = np.median(np.diff(times))
        if tstep == 0:
            thresh_rate = thresh
        else:
            thresh_rate = thresh / np.median(np.diff(times))

        # Set the python time quantity
        time_quantity = pq.second  # Set default
        if "minute" in time_units:
            time_quantity = pq.minute
        elif "hour" in time_units:
            time_quantity = pq.hour
        elif "day" in time_units:
            time_quantity = pq.day

        return thresh_rate

    def get_spike_thresholds(self, values):
        """
        Return the min/max thresholds used for the spike test

        :param values: numpy array of values
        """
        std = np.nanstd(values)
        min_thresh = 1.0 * std
        max_thresh = 2.0 * std
        return min_thresh, max_thresh

    def get_unmasked(self, ncvariable):
        times = self.ncfile.variables["time"][:]
        values = ncvariable[:]
        mask = np.zeros(times.shape[0], dtype=bool)
        if hasattr(values, "mask"):
            mask |= values.mask
        if hasattr(times, "mask"):
            mask |= times.mask
        values = ma.getdata(values[~mask])
        # values = self.normalize_variable(values, ncvariable.units, ncvariable.standard_name)
        return times, values, mask

    def get_values(self, ncvariable):
        times = self.ncfile.variables["time"][:]
        values = ncvariable[:]
        return times, values

    def apply_primary_qc(self, log, ncvariable):
        """
        Applies the primary QC array which is an aggregate of all the other QC
        tests.

        :param netCDF4.Variable ncvariable: NCVariable
        """
        primary_qc_name = "qartod_%s_primary_flag" % ncvariable
        if primary_qc_name not in self.ncfile.variables:
            return
        qcvar = self.ncfile.variables[primary_qc_name]
        qc_variables = self.find_qc_flags(ncvariable)
        vectors = []

        log.info(f"PRIMARY QC FLAG VARS  -- {qc_variables}")
        for qc_variable in qc_variables:
            ncvar = self.ncfile.variables[qc_variable]
            vectors.append(ma.getdata(ncvar[:]))

        flags = qc.qartod_compare(vectors)
        qcvar[:] = flags


def run_qc(log, data, params, ncfile, syntax_test_failed_rows):
    qc = SensorQC(data, params, ncfile, syntax_test_failed_rows)

    for group in geophysical_variables():
        for varname in geophysical_variables()[group]:
            ncvar = geophysical_variables()[group][varname]
            log.info("~" * 75)
            log.info("INSPECTING -- %s | %s", varname, ncvar)
            # log.info("=" * 75)

            if group != "qartod":
                log.info("%s does not need QARTOD", varname)
                continue

            qc.create_variables(log, ncvar, varname)

            for qcvarname in qc.create_qc_variables(log, ncvar, varname):
                log.info("CREATING QC VARIABLE  -- %s", qcvarname)
                qcvar = ncfile.variables[qcvarname]
                log.info("APPLYING QC FOR       -- %s", qcvar.name)
                qc.apply_qc(syntax_test_failed_rows, qcvar)

            qc.apply_primary_qc(log, ncvar)


def main() -> None:
    """Main function"""
    parser = clparser()
    args = parser.parse_args()
    # Arguments
    sensor_file = args.input_file
    num_headerrows = args.header_rows
    param_file = args.param_file
    output_dir = args.output_dir
    log_dir = args.log_dir
    # Options
    verbose = args.verbose
    error_monitoring = args.error_monitoring

    # Validate input arguments
    try:
        num_headerrows = int(num_headerrows)
    except Exception:
        raise ScriptError(ERROR_LEVEL[2], "The [header_rows] argument is invalid.")
    if not os.path.exists(sensor_file):
        raise ScriptError(
            ERROR_LEVEL[2], "The input sensor file does not exist.", sensor_file
        )

    # Abstract some variables based on input arguments
    # Input file basename
    sensor_file_basefn = os.path.basename(sensor_file)
    # NetCDF output
    create_dirs(output_dir)
    ncfile_dir = os.path.join(output_dir, "netcdf")
    create_dirs(ncfile_dir)
    ncfile_basefn = f"{sensor_file_basefn}.nc"
    output_ncfile = os.path.join(ncfile_dir, ncfile_basefn)
    # SQLite database
    db_path = os.path.join(output_dir, "database")
    create_dirs(db_path)
    db_file = os.path.join(db_path, f"{ncfile_basefn}.db")

    # Start logging
    log = start_logging(log_dir)
    if verbose:
        log.addHandler(logging.StreamHandler())
    log_preamble(
        log,
        sensor_file=sensor_file,
        sensor_file_basefn=sensor_file_basefn,
        num_headerrows=num_headerrows,
        parameter_file=param_file,
        output_dir=output_dir,
        ncfile_dir=ncfile_dir,
        ncfile_basefn=ncfile_basefn,
        output_ncfile=output_ncfile,
        db_file=db_file,
        log_dir=log_dir,
        verbose=verbose,
        error_monitoring=error_monitoring,
    )
    log.info("~" * 75)
    log.info("EVALUATING INPUT FILE")
    sensor_obj = FileParser(
        fname=sensor_file, header_rows=num_headerrows, column_delimiter="\t"
    )
    log.info(sensor_obj.filestats())
    log.info("~" * 75)

    # Ingest sensor data as dataframe
    data = MakeDataFrame(sensor_obj)
    data.sensor_dataframe()
    syntax_test_failed_rows = data.syntax_test_failed_rows

    # Ingest QC parameters as dataframe
    param_csv = convert_xlsx(param_file)
    param_obj = FileParser(fname=param_csv, header_rows=0, column_delimiter=",")
    params = MakeDataFrame(param_obj)
    params.parameter_dataframe()
    os.remove(param_csv)

    # Initialize a NetCDF instance with metadata and data dimensions
    ncfile = NetCDF(output_ncfile)
    ncfile.metadata()
    ncfile.operational_variables(data)

    # Run Qartod QC tests
    run_qc(log, data, params, ncfile.rootgrp, syntax_test_failed_rows)

    # --------------------------------
    # Pretty print NetCDF file information, also store sections as variables
    if verbose:
        nc_attrs, nc_dims, nc_vars = ncdump(ncfile.rootgrp, log)
    else:
        nc_attrs, nc_dims, nc_vars = ncdump(ncfile.rootgrp, log, verb=False)

    # --------------------------------

    log.info("~" * 75)
    log.info("Archiving results")
    ds = xr.open_dataset(ncfile.filename)
    df = ds.to_dataframe()
    db_dir = "data/processed/database"
    csv_dir = "data/processed/csv"
    for dir in [db_dir, csv_dir]:
        if not os.path.exists(dir):
            os.mkdir(dir)
    log.info("  SQLite")
    db_path = os.path.join("data/processed/database", os.path.basename(ncfile.filename))
    con = sqlite3.connect(f"{db_path}.db")
    df.to_sql(name="asv_ctd", con=con, if_exists="append")
    csv_path = os.path.join("data/processed/csv", os.path.basename(ncfile.filename))
    log.info("  CSV")
    df.to_csv(f"{csv_path}.csv")

    # --------------------------------

    # Check for errors to report on
    if len(error_list) > 0:
        if verbose:
            print_errors(log, error_list)
        if error_monitoring:
            notifier = Notifier(
                HOST,
                NOTIFIER_EMAIL_FROM,
                NOTIFIER_EMAIL_TO,
                NOTIFIER_EMAIL_CC,
                "Test email",
                print_errors(log, error_list),
            )
            notifier.send()


if __name__ == "__main__":
    main()
