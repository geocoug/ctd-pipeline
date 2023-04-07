# asv-ctd-qc

[![ci](https://github.com/IntegralEnvision/asv-ctd-qc/actions/workflows/ci.yml/badge.svg)](https://github.com/IntegralEnvision/asv-ctd-qc/actions/workflows/ci.yml)

Convert CTD data files from ASCII to NetCDF and flag observations using [QARTOD](https://ioos.noaa.gov/project/qartod/) recommended testing implemented by [ioos_qc](https://github.com/ioos/ioos_qc) :white_check_mark:.

## Setup

### Python (v3.10)

  ```sh
  python -m venv .venv && \  # Create virtual environment
  source ./.venv/bin/activate && \ # Activate the environment
  pip install -r requirements.txt  # Install dependencies
  ```

### Docker

```sh
# Build the Docker image
docker build -t asv-ctd .
```

## Configure QC Test Parameters

IOOS QC test configurations are defined in [config.json](config.json). Below is an example QC configuration for the parameter `pressure`. Keys directly correlate to QC functions and arguments in the [ioos_qc](https://github.com/ioos/ioos_qc) library.

```json
"pressure": {              // Standard CF-safe parameter name
  "axds": {                         // Tests based on the IOOS QC manual
    "valid_range_test": {               // Checks that values are within a min/max range. This is not unlike a `qartod.gross_range_test` with fail and suspect bounds being equal, except that here we specify the inclusive range that should pass instead of the exclusive bounds which should fail
      "valid_span": [-5, 35]                // Values outside the range will FAIL
    }
  },
  "qartod": {                       // Tests based on the IOOS QC manual
    "gross_range_test": {               // Checks that values are within reasonable range bounds.
      "suspect_span": [20, 35],             // Values outside the range will be SUSPECT
       "fail_span": [-5, 35]                // Values outside the range will be FAIL
    },
    "flat_line_test": {                 // Check for consecutively repeated values within a tolerance.
      "tolerance": 0.01,                    // The tolerance that should be exceeded between consecutive values.
      "suspect_threshold": 300,             // The number of seconds within `tolerance` to allow before being flagged as SUSPECT.
      "fail_threshold": 900                 // The number of seconds within `tolerance` to allow before being flagged as FAIL.
    },
    "rate_of_change_test": {            // Checks the first order difference of a series of values to see if there are any values exceeding a threshold defined by the inputs. These are then marked as SUSPECT.
      "threshold": 5                        // A float value representing a rate of change over time, in observation units per second.
    },
    "spike_test": {                     // Check for spikes by checking neighboring data against thresholds.
      "suspect_threshold": 0.01,            // The SUSPECT threshold value, in observations units.
      "fail_threshold": 0.02,               // The SUSPECT threshold value, in observations units.
      "method": "average"                   // ['average'(default),'differential'] optional input to assign the method used to detect spikes.
                                                // "average": Determine if there is a spike at data point n-1 by subtracting the midpoint of n and n-2 and taking the absolute value of this quantity, and checking if it exceeds a low or high threshold.
                                                // "differential": Determine if there is a spike at data point n by calculating the difference between n and n-1 and n+1 and n variation. To considered, (n - n-1)*(n+1 - n) should be smaller than zero (in opposite direction).
    },
    // "attenuated_signal_test": {         // Check for near-flat-line conditions using a range or standard deviation.
    //   "suspect_threshold": 5,               // Any calculated value below this amount will be flagged as SUSPECT. In observations units.
    //   "fail_threshold": 3,                  // Any calculated values below this amount will be flagged as FAIL. In observations units.
    //   "test_period": 15,                    // Length of time to test over in seconds [optional]. Otherwise, will test against entire `inp`.
    //   "min_obs": null,                      // Minimum number of observations in window required to calculate a result [optional]. Otherwise, test will start at beginning of time series. Note: you can specify either `min_obs` or `min_period`, but not both.
    //   "min_period": null,                   // Minimum number of seconds in test_period required to calculate a result [optional]. Otherwise, test will start at beginning of time series. Note: you can specify either `min_obs` or `min_period`, but not both.
    //   "check_type": "range"                 // Either 'std' (default) or 'range', depending on the type of check you wish to perform.
    // },
    // "density_inversion_test": {         // With few exceptions, potential water density will increase with increasing pressure. When vertical profile data is obtained, this test is used to flag as failed T, C, and SP observations, which yield densities that do not sufficiently increase with pressure. A small operator-selected density threshold (DT) allows for micro-turbulent exceptions.
    //   "suspect_threshold": 3,               // A float value representing a maximum potential density(or sigma0) variation to be tolerated, downward density variation exceeding this will be flagged as SUSPECT.
    //   "fail_threshold": 5                   // A float value representing a maximum potential density(or sigma0) variation to be tolerated, downward density variation exceeding this will be flagged as FAIL.
    // },
    "climatology_test": {               // Checks that values are within reasonable range bounds and flags as SUSPECT.
      "zspan": [0, 100]                       // zspan: (optional) Vertical (depth) range, in meters positive down
    }
  }
}
```

## Run QC Checks

See usage notes using the `--help` flag in [asv_ctd_qc.py](./asv_ctd_qc.py).

  ```sh
usage: asv_ctd_qc.py [-h] [-l [LOG]] [-p] [-v] [-c] config input_file output_dir

Evaulate a data file of ASV CTD readings and apply quality assurence checks following QARTOD methods and
assigning data quality flags as appropriate. Transform results into NetCDF format following IC standards.

positional arguments:
  config                Configuration JSON file.
  input_file            Path to the input sensor data file.
  output_dir            Directory for output files.

options:
  -h, --help            show this help message and exit
  -l [LOG], --log [LOG]
                        Path to a log file for script level logging
  -p, --plot            Create an HTML file containing plots of QC flags. Files are stored under a subdirectory of the specified output_dir.
  -v, --verbose         Control the amount of information to display.
  -c, --compliance      Run IOOS compliance checker on compiled NetCDF file.
  ```

### Example QC Run

```sh
docker run -it --rm -v "$(pwd)":/usr/local/app $(docker build -q -t asv-ctd .) python asv_ctd_qc.py -l ./logs/logfile.log -p -c -v config.json ./data/received/2022-10-07T19-45-27.0.txt ./data/processed
```

Output from the above example will produce the following files:

- NetCDF file of sensor data appended with QC flags.
- A compliance check file indicating if any corrective actions are needed for CF compliance.
- PNG files of every variable plotted against each QC check results. PNG files are compiled into a single HTML document for ease of viewing.
- A log file containing script procedures. The same procedures are output to the console using `verbose` mode.

## Plotting QC Flags

Visualize observations and QC flags together using [qc_plots.py](./utils/qc_plots.py). It is also built in to [asv_ctd_qc.py](./asv_ctd_qc.py) as an optional argument.

```sh
usage: qc_plots.py [-h] [-v] ncfile outdir

Create plots of ASV CTD cast data with QARTOD flags.

positional arguments:
  ncfile         Path to the input NetCDF file.
  outdir         Path to save output files.

options:
  -h, --help     show this help message and exit
  -v, --verbose  Control the amount of information to display.
```

## Dumping NetCDF Contents

The `netcdf` package shows details of file contents, which can be downloaded using the `homebrew` or `apt` package managers. Likewise, use the [ncdump.py](./utils/ncdump.py) utility with Python. [ncdump.py](./utils/ncdump.py) is built in to [asv_ctd_qc.py](./asv_ctd_qc.py) as an optional argument.

## CF Compliance Checks

> The [IOOS Compliance Checker](https://github.com/ioos/compliance-checker) is a python based tool for data providers to check for completeness and community standard compliance of local or remote netCDF files against CF and ACDD file standards. The python module can be used as a command-line tool or as a library that can be integrated into other software.

Compliance checks can be run using the command line utility, or using Python via [compliance.py](./utils/compliance.py). Compliance checks are also built in to [asv_ctd_qc.py](./asv-ctd_qc.py) as an optional argument.

## Notes

### ioos_qc

- [ioos_qc](https://github.com/ioos/ioos_qc) was cloned to the current repository on version 2.1.0 and modified for project requirements. The following revisions were made:

#### qartod.py

1. [time_interval](https://github.com/IntegralEnvision/asv-ctd-qc/commit/a249dd4ee84f719696fb31ecd6eabd9edd0f6a33#diff-32c09032f00f303300ace35369debee33af51ceb355defcce878c489bdc3af6aR646) calculation. CTD collection times can be < 1 second apart causing the number of observations to appear as 0.

#### streams.py

1. Updated `self.lat_column` from `lat` to `latitude` and `self.lon_column` from `lon` to `longitude` in [`PandasStream().__init__()`](https://github.com/ioos/ioos_qc/blob/093935e0f2c21a6a585bda5a194fc7a2c7aedd76/ioos_qc/streams.py#L49)

### Souce Data String Configuration (tab delimited)

- Data / Time
- Turbidity (ftu - formazin turbidity unit)
- Pressure (dbar - decibar)
- Temperature (°C)
- Dissolved Oxygen (% saturation)
- Altitude (meters)
- Conductivity (mS/cm - millisiemens per centimeter)
- Fluorometer (ug/L - microgram per liter)
- pH
- Calculated Salinity (psu - practical salinity unit)
- Calculated Depth (meters)
- Latitude (degrees)
- Longitude (degrees)

## Help

- [QARTOD](https://ioos.noaa.gov/project/qartod/)
- [Intro to NetCDF](https://adyork.github.io/python-oceanography-lesson/17-Intro-NetCDF/index.html)
- [Writing NetCDF](https://www.earthinversion.com/utilities/Writing-NetCDF4-Data-using-Python/)
- [Dimensions](http://www.bic.mni.mcgill.ca/users/sean/Docs/netcdf/guide.txn_12.html)
- [Western Indian Ocean (WIO) Workshop](https://github.com/MathewBiddle/WIO_workshop)

## Abstract

**Data Management Workflow for Autonomous Surface Vehicle CTD Data: A Python-based Approach for Near Real-time Processing and Quality Control**

The increasing deployment of autonomous surface vehicles (ASVs) for oceanographic research has led to a growing need for efficient and reliable data management workflows. This paper presents a Python-based workflow for handling and processing CTD (Conductivity, Temperature, and Depth) data collected by ASVs in near real-time. The workflow focuses on the conversion of raw CTD data files from ASCII format to the more versatile NetCDF format, which facilitates data sharing and interoperability among researchers and institutions. Additionally, the workflow incorporates quality control measures based on the Quality Assurance of Real-Time Oceanographic Data (QARTOD) recommendations provided by the Integrated Ocean Observing System (IOOS). These measures ensure that the processed CTD data is accurate and reliable for further analysis and interpretation. The proposed workflow demonstrates a robust and efficient approach to managing ASV CTD data, streamlining the data processing pipeline and promoting the use of high-quality oceanographic data in research and decision-making processes.
