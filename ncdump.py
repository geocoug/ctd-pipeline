import logging

import xarray as xr


def ncdump(nc_fid: xr.Dataset, log: logging.Logger, verb=True) -> tuple:
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
            log.info("\t\ttype: %s" % repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                log.info(
                    "\t\t%s: %s"
                    % (ncattr, repr(nc_fid.variables[key].getncattr(ncattr))),
                )
        except KeyError:
            log.info("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_fid = nc_fid.rootgrp
    nc_attrs = nc_fid.ncattrs()
    if verb:
        log.info("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            log.info(f"\t{nc_attr}: {repr(nc_fid.getncattr(nc_attr))}")
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
