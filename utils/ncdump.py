#!/usr/bin/env python

import argparse
import logging

import netCDF4 as nc  # noqa N813
import xarray as xr

logger = logging.getLogger("__main__")


class NetCDF:
    """Create a NetCDF object and provide methods for reading
    self-describing attributes.
    """

    def __init__(
        self: "NetCDF",
        filename: str,
        data_model: str = "NETCDF4",
    ) -> None:
        """Create a NetCDF object and provide methods for reading
        self-describing attributes.

        Args:
        ----
            ncfile (str): NetCDF file.
            data_model (str, optional): NC data model. Defaults to "NETCDF4".
        """
        self.data_model = data_model
        self.filename = filename
        self.rootgrp = nc.Dataset(self.filename, "r", format=self.data_model)


def ncdump(
    nc_fid: xr.Dataset,
) -> tuple:
    """Ncdump outputs dimensions, variables and their attribute information.

    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object.

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
        """Prints the NetCDF file attributes for a given key.

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

    logger.info("LOGGING OUTPUT NETCDF ATTRIBUTES")

    # NetCDF global attributes
    nc_fid = nc_fid.rootgrp
    nc_attrs = nc_fid.ncattrs()
    logger.info("NetCDF Global Attributes:")
    for nc_attr in nc_attrs:
        logger.info(f"\t{nc_attr}: {repr(nc_fid.getncattr(nc_attr))}")
    nc_dims = list(nc_fid.dimensions)  # list of nc dimensions
    # Dimension shape information.
    logger.info("NetCDF dimension information:")
    for dim in nc_dims:
        logger.info("\tName: %s" % dim)
        logger.info("\t\tsize: %s" % len(nc_fid.dimensions[dim]))
        print_ncattr(dim)
    # Variable information.
    nc_vars = list(nc_fid.variables)  # list of nc variables
    logger.info("NetCDF variable information:")
    for var in nc_vars:
        if var not in nc_dims:
            logger.info("\tName: %s" % var)
            logger.info("\t\tdimensions: %s" % nc_fid.variables[var].dimensions)
            logger.info("\t\tsize: %s" % nc_fid.variables[var].size)
            print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars


def clparser() -> argparse.ArgumentParser:
    """Create a parser to handle input arguments."""
    desc_msg = """Dump contents of a NetCDF file."""
    parser = argparse.ArgumentParser(description=desc_msg)
    parser.add_argument("ncfile", help="Path to the input NetCDF file.")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        help="Control the amount of information to display.",
    )
    return parser


def main() -> None:
    """Main entrypoint."""
    parser = clparser()
    args = parser.parse_args()
    ncfile = args.ncfile
    verbose = args.verbose
    if verbose:
        logger.addHandler(logging.StreamHandler())
        logger.setLevel(logging.INFO)
    ncdump(nc_fid=NetCDF(ncfile))


if __name__ == "__main__":
    main()
