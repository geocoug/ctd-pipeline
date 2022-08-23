#!/usr/bin/env python
# coding=utf-8

import datetime
import logging
import os
import sys

import xarray as xr

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


def ncdump(nc_fid: xr.Dataset, verb=True) -> tuple:
    """
    Ncdump outputs dimensions, variables and their attribute information.

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

    def print_ncattr(key: str) -> None:
        """
        Prints the NetCDF file attributes for a given key.

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            logger.info("\t\ttype: %s" % repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                logger.info(
                    "\t\t%s: %s"
                    % (ncattr, repr(nc_fid.variables[key].getncattr(ncattr))),
                )
        except KeyError:
            logger.info("\t\tWARNING: %s does not contain variable attributes" % key)

    if verb:
        logger.addHandler(logging.StreamHandler())
    # NetCDF global attributes
    nc_fid = nc_fid.rootgrp
    nc_attrs = nc_fid.ncattrs()
    if verb:
        logger.info("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            logger.info(f"\t{nc_attr}: {repr(nc_fid.getncattr(nc_attr))}")
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        logger.info("NetCDF dimension information:")
        for dim in nc_dims:
            logger.info("\tName: %s" % dim)
            logger.info("\t\tsize: %s" % len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        logger.info("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                logger.info("\tName: %s" % var)
                logger.info("\t\tdimensions: %s" % nc_fid.variables[var].dimensions)
                logger.info("\t\tsize: %s" % nc_fid.variables[var].size)
                print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars
