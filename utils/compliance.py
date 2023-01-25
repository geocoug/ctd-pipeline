#!/usr/bin/env python
import os

from compliance_checker.runner import CheckSuite, ComplianceChecker


class Compliance:
    """Shorthand for ComplianceChecker"""

    def __init__(self: "Compliance", convention: str) -> None:
        """Initialize a compliance checker using a single convention.

        Args:
            convention (str): CF convention.
            See the [docs](https://github.com/ioos/compliance-checker) for more info.
        """
        self.convention = convention
        self.check_suite = CheckSuite()
        self.check_suite.load_all_available_checkers()

    def run_checker(
        self: "Compliance",
        ncfile: str,
        output_filename: str,
        verbose: bool = False,
    ) -> None:
        """Shorthand to run ComplianceChecker

        Args:
            ncfile (str): NetCDF file to evaluate for compliance.
            output_filename (str): Compliance output filename.
            verbose (bool, optional): Verbose mode. Defaults to False.
        """
        ComplianceChecker.run_checker(
            ds_loc=ncfile,
            checker_names=[self.convention.replace("-", ":").lower()],
            verbose=verbose,
            criteria="normal",
            output_filename=output_filename,
            output_format=os.path.splitext(output_filename)[-1].replace(".", ""),
        )
