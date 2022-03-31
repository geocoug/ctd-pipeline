#! /usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import csv
import datetime
import logging
import os

import netCDF4 as nc
import numpy as np
import numpy.ma as ma
import pandas as pd
import quantities as pq
import xlsx2csv as xl

from ioos_qc import argo, axds
from ioos_qc import qartod as qc
from ioos_qc import utils

logdir = "logs"
if not os.path.exists(logdir):
    os.mkdir(logdir)
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s : %(msecs)04d : %(levelname)s : %(message)s",
                    filename=(os.path.join(logdir, 'log_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + ".log")),
                    filemode='w',
                    datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger(__name__)

def ncdump(nc_fid, verb=True):
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
                log.info('\t\t%s: %s' % (
                    ncattr, repr(nc_fid.variables[key].getncattr(ncattr))))
        except KeyError:
            log.info("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        log.info("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            log.info('\t%s: %s' % (nc_attr, repr(nc_fid.getncattr(nc_attr))))
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
                log.info('\tName: %s' % var)
                log.info("\t\tdimensions: %s" % nc_fid.variables[var].dimensions)
                log.info("\t\tsize: %s" % nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars


def geophysical_variables():
    '''
    Returns a dict of variable names matching the geophysical variable's
    standard names
    '''
    standard_name_lookup = {
        'qartod': {
            'TEMPERATURE;C': 'sea_water_temperature',
            'CONDUCTIVITY;MS/CM': 'sea_water_electrical_conductivity',
            'PRESSURE;DBAR': 'sea_water_pressure',
            'Calc. SALINITY; PSU': 'sea_water_practical_salinity',
            'PH;PH': 'sea_water_ph',
            'DISSOLVED OXYGEN;SAT%': 'sea_water_dissolved_oxygen',
            'FLUOROMETER (C);UG/L': 'sea_water_fluorescence',
            'TURBIDITY;FTU': 'sea_water_turbidity',
            'ALTITUDE;M': 'altitude',
            'Calc. DEPTH;M': 'depth'
        },
        'operational': {
            'Date / Time': 'datetime',
            # 'ALTITUDE;M'            : 'altitude',
            # 'Calc. DEPTH;M'         : 'depth',
            'LATITUDE;DEG': 'latitude',
            'LONGITUDE;DEG': 'longitude'
        }
    }
    return standard_name_lookup


file_error_list = []


class FileError(Exception):
    """Handle compilation of file manipulation errors into a
    somewhat readable format for debugging purposes."""
    def __init__(self, errmsg, infile=None):
        self.errmsg = errmsg
        self.infile = infile
        file_error_list.append((self.errmsg, self.infile))


def clparser():
    """Create a parser to handle input arguments and displaying
    a script specific help message."""
    desc_msg = """Evaulate a data file of ASV CTD readings
        and applies quality assurence
        checks following QARTOD methods and assigning data
        quality flags as appropriate. Transform
        results into NetCDF format following IC standards."""
    parser = argparse.ArgumentParser(description=desc_msg)
    parser.add_argument('input_file',
                        help="Path to the input sensor data file.")
    parser.add_argument('header_rows',
                        help="""Number of rows preceeding the row
                        containing column headers.""")
    parser.add_argument('sensor_id',
                        help='Identification tag of the sensor.')
    parser.add_argument('sensor_sn',
                        help='Serial number of the sensor.')
    parser.add_argument('output_file',
                        help='Path for the output NetCDF file.')
    parser.add_argument('param_file',
                        help='Path to sensor threshold parameters file.')
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False,
                        help="Control the amount of information to display.")
    return parser


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
        return self


class FileParser:
    """Parse input file attributes"""
    def __init__(self, fname, header_rows, column_delimiter):
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
            FileError("Can't open file to determine format.", self.filename)
        try:
            if self.encoding:
                reader = UnicodeReader(open(self.filename, "rt"), dialect, self.encoding)
            else:
                reader = csv.reader(open(self.filename, "rt"), dialect)
        except Exception:
            FileError("Can't open file to read data.", self.filename)
        if self.headrows > 0:
            try:
                for row in range(self.headrows):
                    next(reader)
            except Exception:
                FileError("Can't read column header line.", self.filename)
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
        fstats = "Filename: %s\nFile rows: %d\nHeader rows: %d\nColumn row: %d\nData rows: %d\n" % (self.filename, self.filerows, self.headrows, self.colrow, self.datarows)
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
        for (idx, row) in enumerate(self.reader):
            row = list(filter(None, row))
            if len(row) < len(self.cols):
                # MISSING VALUES SHOULD BE AUTOMATIC FAIL FLAG. DO NOT REPLACE WITH 0
                FileError("The input file is missing columns on line %s. Missing values will be replaced with 0." % (idx + self.fobj.colrow + 1), self.fobj.filename)
                for c in range(len(self.cols) - len(row)):
                    row.append('0.0')
            if len(row) > len(self.cols):
                FileError("The input file is missing column headers on line %s. Columns with no headers will be removed." % (idx + self.fobj.colrow + 1), self.fobj.filename)
                row = row[:-(len(row) - len(self.cols))]
            self.data_array.append(row)
        self.df = pd.DataFrame(self.data_array, columns=self.cols)
        """Set data types of the DataFrame"""
        for col in range(len(self.cols)):
            if col == 0:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype('datetime64')
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]] - datetime.datetime(2010, 1, 1, 0, 0)
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].dt.total_seconds()
            else:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype('float')
                
    def parameter_dataframe(self):
        """Create a Pandas DataFrame from an array
        of arrays retrieved from the reader"""
        self.data_array = []
        for (idx, row) in enumerate(self.reader):
            self.data_array.append(list(filter(None, row)))
        self.df = pd.DataFrame(self.data_array, columns=self.cols)
        """Set data types of the DataFrame"""
        for col in range(len(self.cols)):
            if 'value' in self.cols[col]:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype('float')
            else:
                self.df[self.df.columns[col]] = self.df[self.df.columns[col]].astype('str')

def operational_test():
    pass


class SensorQC(object):
    def __init__(self, data, config_params, ncfile):
        self.data = data
        self.params = config_params
        self.ncfile = ncfile

    def create_variables(self, ncvariable, varname):
        """
        Returns a list of variable names for the newly created variables
        """
        name = ncvariable
        standard_name = ncvariable
        dims = 'time'
        units = varname.split(';')[1].strip()
        source_name = varname
        log.info("Creating variables for %s", name)
        variables = []

        variable_name = name
        ncvar = self.ncfile.createVariable(variable_name, np.float64, dims, fill_value=np.int8(9))
        ncvar.source_name = source_name
        ncvar.units = units
        ncvar.standard_name = standard_name
        ncvar.long_name = standard_name
        ncvar[:] = self.data.df[source_name]
        variables.append(variable_name)
        return variables


    def find_qc_flags(self, ncvariable):
        '''
        Returns a list of non-GliderDAC QC flags associated with a variable

        :param netCDF4.Variable ncvariable: Variable to get the status flag
                                            variables for
        '''
        valid_variables = []
        ancillary_variables = getattr(ncvariable, 'ancillary_variables', None)
        if isinstance(ancillary_variables, str):
            ancillary_variables = ancillary_variables.split(' ')
        else:
            return []
        for varname in ancillary_variables:
            if varname not in self.ncfile.variables:
                log.warning("%s defined as ancillary variable but doesn't exist", varname)
                continue
            if varname.endswith('_qc'):
                valid_variables.append(varname)
            if 'status_flag' in getattr(self.ncfile.variables[varname], 'standard_name', ''):
                valid_variables.append(varname)

        return valid_variables
    
    
    def create_qc_variables(self, ncvariable, varname):
        """
        Returns a list of variable names for the newly created variables for QC flags
        """
        name = ncvariable
        standard_name = ncvariable
        dims = 'time'
        log.info("Creating QARTOD variables for %s", name)

        templates = {
            'gap': {
                'name': 'qartod_%(name)s_gap_flag',
                'long_name': 'QARTOD Gap Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'gap',
                'dac_comment': 'ioos_qartod'
            },
            'valid_range': {
                'name': 'qartod_%(name)s_valid_range_flag',
                'long_name': 'QARTOD Valid Range Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'valid_range',
                'dac_comment': 'ioos_qartod'
            },
            'location': {
                'name': 'qartod_%(name)s_location_flag',
                'long_name': 'QARTOD Location Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'location',
                'dac_comment': 'ioos_qartod'
            },
            # 'syntax': {
            #     'name': 'qartod_%(name)s_syntax_flag',
            #     'long_name': 'QARTOD Syntax Test for %(standard_name)s',
            #     'standard_name': '%(standard_name)s status_flag',
            #     'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
            #     'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
            #     'qartod_test': 'syntax',
            #     'dac_comment': 'ioos_qartod'
            # },
            # 'climatological': {
            #     'name': 'qartod_%(name)s_climatological_flag',
            #     'long_name': 'QARTOD Climatological Test for %(standard_name)s',
            #     'standard_name': '%(standard_name)s status_flag',
            #     'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
            #     'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
            #     'qartod_test': 'climatological',
            #     'dac_comment': 'ioos_qartod'
            # },
            'flat_line': {
                'name': 'qartod_%(name)s_flat_line_flag',
                'long_name': 'QARTOD Flat Line Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'flat_line',
                'dac_comment': 'ioos_qartod'
            },
            'gross_range': {
                'name': 'qartod_%(name)s_gross_range_flag',
                'long_name': 'QARTOD Gross Range Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'gross_range',
                'dac_comment': 'ioos_qartod'
            },
            'rate_of_change': {
                'name': 'qartod_%(name)s_rate_of_change_flag',
                'long_name': 'QARTOD Rate of Change Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'rate_of_change',
                'dac_comment': 'ioos_qartod'
            },
            'spike': {
                'name': 'qartod_%(name)s_spike_flag',
                'long_name': 'QARTOD Spike Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'spike',
                'dac_comment': 'ioos_qartod'
            },
            'pressure': {
                'name': 'qartod_%(name)s_pressure_flag',
                'long_name': 'QARTOD Pressure Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'pressure',
                'dac_comment': 'ioos_qartod'
            },
            'attenuated_signal': {
                'name': 'qartod_%(name)s_attenuated_signal_flag',
                'long_name': 'QARTOD Attenuated Signal Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'attenuated_signal',
                'dac_comment': 'ioos_qartod'
            },
            'density_inversion': {
                'name': 'qartod_%(name)s_density_inversion_flag',
                'long_name': 'QARTOD Density Inversion Test for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'qartod_test': 'density_inversion',
                'dac_comment': 'ioos_qartod'
            },
            'primary': {
                'name': 'qartod_%(name)s_primary_flag',
                'long_name': 'QARTOD Primary Flag for %(standard_name)s',
                'standard_name': '%(standard_name)s status_flag',
                'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
                'flag_meanings': 'GOOD NOT_EVALUATED SUSPECT BAD MISSING',
                'dac_comment': 'ioos_qartod_primary'
            }
        }

        qcvariables = []

        for tname, template in list(templates.items()):
            if tname == 'pressure' and standard_name != 'sea_water_pressure':
                continue
            variable_name = template['name'] % {'name': name}

            if variable_name not in geophysical_variables()['qartod']:
                ncvar = self.ncfile.createVariable(variable_name, np.int8, dims, fill_value=np.int8(9))
            else:
                ncvar = self.ncfile.variables[variable_name]

            ncvar.units = '1'
            ncvar.standard_name = template['standard_name'] % {'standard_name': standard_name}
            ncvar.long_name = template['long_name'] % {'standard_name': standard_name}
            ncvar.flag_values = template['flag_values']
            ncvar.flag_meanings = template['flag_meanings']
            # ncvar.references = template['references']
            ncvar.dac_comment = template['dac_comment']
            if 'qartod_test' in template:
                ncvar.qartod_test = template['qartod_test']
            qcvariables.append(variable_name)

        return qcvariables


    def apply_qc(self, ncvariable):
        """
        Applies QC to a qartod variable
        """
        qc_tests = {
            'gap': utils.check_timestamps,
            'valid_range': axds.valid_range_test,
            'location': qc.location_test,
            # 'syntax': ??,
            # 'climatological': qc.climatology_test,
            'flat_line': qc.flat_line_test,
            'gross_range': qc.gross_range_test,
            'rate_of_change': qc.rate_of_change_test,
            'spike': qc.spike_test,
            'attenuated_signal': qc.attenuated_signal_test,
            'pressure': argo.pressure_increasing_test,
            'density_inversion': qc.density_inversion_test
        }

        qartod_test = getattr(ncvariable, 'qartod_test', None)
        if not qartod_test:
            return
        standard_name = getattr(ncvariable, 'standard_name').split(' ')[0]
        parent = self.ncfile.get_variables_by_attributes(standard_name=standard_name)[0]
        times, values, mask = self.get_unmasked(parent)
        # There's no data to QC
        if len(values) == 0:
            return

        # If the config isn't set for this test, don't run it.
        # if qartod_test not in self.config[standard_name]:
            # return

        qa_config = self.params.df[self.params.df['sensor'] == getattr(parent, 'source_name')]
        op_config = self.params.df[self.params.df['sensor'] == 'OPERATOR']

        test_params = {}

        if qartod_test == 'gap':
            test_params['times'] = ma.getdata(times[~mask])

        if qartod_test == 'location':
            test_params['lat'] = self.ncfile.get_variables_by_attributes(standard_name='latitude')[0][:]
            test_params['lon'] = self.ncfile.get_variables_by_attributes(standard_name='longitude')[0][:]
            test_params['bbox'] = tuple([       
                op_config.loc[op_config['parameter'] == 'long_min', 'parameter_value'].iloc[0],
                op_config.loc[op_config['parameter'] == 'lat_min', 'parameter_value'].iloc[0],
                op_config.loc[op_config['parameter'] == 'long_max', 'parameter_value'].iloc[0],
                op_config.loc[op_config['parameter'] == 'lat_max', 'parameter_value'].iloc[0]
            ])

        if qartod_test in 'rate_of_change':
            test_params['inp'] = values
            time_units = self.ncfile.variables['time'].units
            dates = np.array(nc.num2date(times, time_units), dtype='datetime64[ms]')
            test_params['tinp'] = dates
            n_dev = qa_config.loc[qa_config['parameter'] == 'n_deviation', 'parameter_value'].iloc[0]
            test_params['threshold'] = self.get_rate_of_change_threshold(values, times, time_units, n_dev)

        if qartod_test == 'spike':
            test_params['inp'] = values
            test_params['suspect_threshold'] = qa_config.loc[qa_config['parameter'] == 'spike_low', 'parameter_value'].iloc[0]
            test_params['fail_threshold'] = qa_config.loc[qa_config['parameter'] == 'spike_high', 'parameter_value'].iloc[0]

        if qartod_test == 'gross_range':
            test_params['inp'] = values
            test_params['fail_span'] = tuple([
                qa_config.loc[qa_config['parameter'] == 'sensor_min', 'parameter_value'].iloc[0],
                qa_config.loc[qa_config['parameter'] == 'sensor_max', 'parameter_value'].iloc[0]
            ])
            test_params['suspect_span'] = tuple([
                qa_config.loc[qa_config['parameter'] == 'climate_min', 'parameter_value'].iloc[0],
                qa_config.loc[qa_config['parameter'] == 'climate_max', 'parameter_value'].iloc[0]
            ])
            
        if qartod_test == 'valid_range':
            test_params['inp'] = values
            test_params['valid_span'] = tuple([
                qa_config.loc[qa_config['parameter'] == 'sensor_min', 'parameter_value'].iloc[0],
                qa_config.loc[qa_config['parameter'] == 'sensor_max', 'parameter_value'].iloc[0]
            ])
            
        # if qartod_test == 'climatological':
        #     test_params['low_thresh'] = qa_config.loc[qa_config['parameter'] == 'climate_min', 'parameter_value'].iloc[0]
        #     test_params['high_thresh'] = qa_config.loc[qa_config['parameter'] == 'climate_max', 'parameter_value'].iloc[0]

        if qartod_test == 'flat_line':
            test_params['inp'] = values
            test_params['tinp'] = ma.getdata(times[~mask])
            test_params['suspect_threshold'] = qa_config.loc[qa_config['parameter'] == 'rep_cnt_susp', 'parameter_value'].iloc[0]
            test_params['fail_threshold'] = qa_config.loc[qa_config['parameter'] == 'rep_cnt_fail', 'parameter_value'].iloc[0]
            test_params['tolerance'] = qa_config.loc[qa_config['parameter'] == 'eps', 'parameter_value'].iloc[0]
            
        if qartod_test == 'pressure':
            test_params['inp'] = values
            
        if qartod_test == 'density_inversion':
            test_params['inp'] = values
            test_params['zinp'] = self.data.df['PRESSURE;DBAR'].to_numpy()
            
        if qartod_test == 'attenuated_signal':
            test_params['inp'] = values
            test_params['tinp'] = ma.getdata(times[~mask])
            test_params['suspect_threshold'] = qa_config.loc[qa_config['parameter'] == 'climate_min', 'parameter_value'].iloc[0]
            test_params['fail_threshold'] = qa_config.loc[qa_config['parameter'] == 'climate_max', 'parameter_value'].iloc[0]
            test_params['min_obs'] = qa_config.loc[qa_config['parameter'] == 'rep_cnt_susp', 'parameter_value'].iloc[0]

        qc_flags = qc_tests[qartod_test](**test_params)
        ncvariable[~mask] = qc_flags

        for test_param in test_params:
            if test_param in ('inp', 'tinp', 'zinp'):
                continue
            ncvariable.setncattr(test_param, test_params[test_param])

    def get_rate_of_change_threshold(self, values, times,
                                     time_units='seconds since 1970-01-01T00:00:00Z',
                                     n_dev=3):
        """
        Return the threshold used for the rate of change test

        :param values: numpy array of values
        :param times: numpy array of times
        :param time_units: string defining time units
        """
        n_dev = n_dev   # Set to 3 standard deviations
        std = np.nanstd(values)
        thresh = n_dev * std
        tstep = np.median(np.diff(times))
        if tstep == 0:
            thresh_rate = thresh
        else:
            thresh_rate = thresh / np.median(np.diff(times))

        # Set the python time quantity
        time_quantity = pq.second  # Set default
        if 'minute' in time_units:
            time_quantity = pq.minute
        elif 'hour' in time_units:
            time_quantity = pq.hour
        elif 'day' in time_units:
            time_quantity = pq.day
            
        return thresh_rate

    def get_spike_thresholds(self, values):
        '''
        Return the min/max thresholds used for the spike test

        :param values: numpy array of values
        '''
        std = np.nanstd(values)
        min_thresh = 1.0 * std
        max_thresh = 2.0 * std
        return min_thresh, max_thresh

    def get_unmasked(self, ncvariable):
        times = self.ncfile.variables['time'][:]
        values = ncvariable[:]
        mask = np.zeros(times.shape[0], dtype=bool)
        if hasattr(values, 'mask'):
            mask |= values.mask
        if hasattr(times, 'mask'):
            mask |= times.mask
        values = ma.getdata(values[~mask])
        # values = self.normalize_variable(values, ncvariable.units, ncvariable.standard_name)
        return times, values, mask

    def get_values(self, ncvariable):
        times = self.ncfile.variables['time'][:]
        values = ncvariable[:]
        return times, values

    def apply_primary_qc(self, ncvariable):
        '''
        Applies the primary QC array which is an aggregate of all the other QC
        tests.

        :param netCDF4.Variable ncvariable: NCVariable
        '''
        primary_qc_name = 'qartod_%s_primary_flag' % ncvariable
        if primary_qc_name not in self.ncfile.variables:
            return

        qcvar = self.ncfile.variables[primary_qc_name]
        # Only perform primary QC on variables created by DAC
        if getattr(qcvar, 'dac_comment', '') != 'ioos_qartod_primary':
            return

        qc_variables = self.find_qc_flags(ncvariable)
        vectors = []

        for qc_variable in qc_variables:
            ncvar = self.ncfile.variables[qc_variable]
            if getattr(ncvar, 'dac_comment', '') != 'ioos_qartod':
                continue
            log.info("Using %s in primary QC for %s", qc_variable, primary_qc_name)
            vectors.append(ma.getdata(ncvar[:]))

        log.info("Applying QC for %s", primary_qc_name)
        flags = qc.qc_compare(vectors)
        qcvar[:] = flags


class NetCDF:
    """Create a NetCDF object and provide methods for adding self-describing attributes"""
    def __init__(self, outfile, data_model="NETCDF4"):
        self.filename = outfile
        self.data_model = data_model
        self.rootgrp = nc.Dataset(self.filename, "w", format=self.data_model)

    def metadata(self, sensor_id, sensor_sn):
        self.rootgrp.title = "Hypoxia Monitoring Data"
        self.rootgrp.description = "NOAA ASV sensor parameter QA"
        self.rootgrp.history = "Created %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.rootgrp.source = "L3"
        self.rootgrp.conventions = "CF-1.6"
        self.rootgrp.author = "Integral Consulting Inc."
        self.rootgrp.sensor_id = sensor_id
        self.rootgrp.sensor_sn = sensor_sn

    def operational_variables(self, data):
        tim_dim = self.rootgrp.createDimension("time", None)
        tim = self.rootgrp.createVariable("time", np.dtype("float64").char, ("time",))
        tim.long_name = "Time (PST) as Seconds Since 2010-01-01 00:00:00"
        tim.standard_name = "time"
        tim.units = "seconds since 2010-01-01 00:00:00"
        tim.time_zone = "PST"
        tim.calendar = "standard"
        tim[:] = data.df['Date / Time']

        lat = self.rootgrp.createVariable("lat", np.dtype("float64").char, ("time",))
        lat.long_name = "Latitude degrees"
        lat.standard_name = "latitude"
        lat.units = "degrees_north"
        lat[:] = data.df['LATITUDE;DEG']

        lon = self.rootgrp.createVariable("lon", np.dtype("float64").char, ("time",))
        lon.long_name = "Longitude degrees"
        lon.standard_name = "longitude"
        lon.units = "degrees_east"
        lon[:] = data.df['LONGITUDE;DEG']

        # alt = self.rootgrp.createVariable("altitude", np.dtype("float64").char, ("time",))
        # alt.long_name = "Altitude..."
        # alt.standard_name = "altitude"
        # alt.units = "meters"
        # alt[:] = data.df['ALTITUDE;M']

        # depth = self.rootgrp.createVariable("depth", np.dtype("float64").char, ("time",))
        # depth.long_name = "Calculated depth below water surface"
        # depth.standard_name = "depth"
        # depth.units = "meters"
        # depth[:] = data.df['Calc. DEPTH;M']

        # temp = self.rootgrp.createVariable("sea_water_temperature", np.dtype("float64").char, ("time",))
        # temp.long_name = "Temperature"
        # temp.standard_name = "sea_water_temperature"
        # temp.units = "celcius"
        # temp[:] = data.df['TEMPERATURE;C']


def run_qc(data, params, ncfile):
    qc = SensorQC(data, params, ncfile)

    for group in geophysical_variables():
        for varname in geophysical_variables()[group]:
            ncvar = geophysical_variables()[group][varname]
            log.info("~" * 75)
            log.info("INSPECTING -- %s | %s", varname, ncvar)
            log.info("=" * 75)

            if group != 'qartod':
                log.info("%s does not need QARTOD", varname)
                continue

            qc.create_variables(ncvar, varname)

            for qcvarname in qc.create_qc_variables(ncvar, varname):
                log.info("CREATING QC VARIABLE  -- %s", qcvarname)
                qcvar = ncfile.variables[qcvarname]
                log.info("APPLYING QC FOR       -- %s", qcvar.name)
                qc.apply_qc(qcvar)

        # qc.apply_primary_qc(ncvar)


def log_preamble(**kwargs):
    """Write a preamble to a log file consisting of key:value pairs."""
    log.info("start time: %s %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"), datetime.datetime.now().astimezone().tzinfo))
    if len(kwargs) > 0:
        for arg in kwargs:
            log.info("%s: %s" % (arg, kwargs[arg]))


def convert_xlsx(file):
    """Convert an XLSX file to CSV."""
    if not os.path.exists(file):
        raise FileError("The input file does not exist.", file)
    filename_noext = os.path.splitext(os.path.basename(os.path.abspath(file)))[0]
    output_csv = os.path.join(os.path.dirname(os.path.abspath(file)),
                              f"{filename_noext}.csv")
    xl.Xlsx2csv(file).convert(output_csv, sheetid=2)
    return output_csv


def main():
    """Main function"""
    parser = clparser()
    args = parser.parse_args()
    # Input arguments
    sensor_file = args.input_file
    num_headerrows = args.header_rows
    sensor_name = args.sensor_id
    sensor_sn = args.sensor_sn
    netcdf_output = args.output_file
    param_file = args.param_file
    # Input options    
    verbose = args.verbose
    if verbose:
        log.addHandler(logging.StreamHandler())

    log_preamble(
        sensor_file=sensor_file,
        num_headerrows=num_headerrows,
        sensor_name=sensor_name,
        sensor_sn=sensor_sn,
        netcdf_output=netcdf_output,
        parameter_file=param_file)

    try:
        num_headerrows = int(num_headerrows)
    except Exception:
        raise FileError("The {Header Rows} argument is invalid.")
    if not os.path.exists(sensor_file):
        raise FileError("The input sensor file does not exist.", sensor_file)

    sensor_obj = FileParser(sensor_file, num_headerrows, column_delimiter="\t")
    data = MakeDataFrame(sensor_obj)
    data.sensor_dataframe()

    param_csv = convert_xlsx(param_file)

    param_obj = FileParser(fname=param_csv, header_rows=0, column_delimiter=",")
    params = MakeDataFrame(param_obj)
    params.parameter_dataframe()
    os.remove(param_csv)

    ncfile = NetCDF(netcdf_output)

    ncfile.metadata(sensor_name, sensor_sn)

    ncfile.operational_variables(data)

    run_qc(data, params, ncfile.rootgrp)


    # print(ncfile.rootgrp.dimensions)
    # print(ncfile.rootgrp.variables.keys())
    # print(ncfile.rootgrp.variables['altitude'].shape)
    # print(ncfile.rootgrp['sea_water_dissolved_oxygen'].dtype)
    # print(ncfile.rootgrp['sea_water_dissolved_oxygen'].ndim)
    # print(ncfile.rootgrp)
    # print(ncfile.rootgrp.variables)
    # print(ncfile.rootgrp.variables['sea_water_temperature'][:])
    # print(ncfile.rootgrp.variables['time'][:][5])
    # print(ncfile.rootgrp.variables['qartod_sea_water_temperature_syntax_flag'][:])

    # --------------------------------
    # Pretty print NetCDF file information, also store sections as variables
    if verbose:
        nc_attrs, nc_dims, nc_vars = ncdump(ncfile.rootgrp)
    else:
        nc_attrs, nc_dims, nc_vars = ncdump(ncfile.rootgrp, verb=False)
    # --------------------------------

    if len(file_error_list) != 0:
        print('')
        for err in file_error_list:
            print(err)

    log.info("=" * 75)
    log.info("Completion time: %s %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S %p"), datetime.datetime.now().astimezone().tzinfo))


if __name__ == "__main__":
    main()
