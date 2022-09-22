# asv-ctd-qa

[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/IntegralEnvision/asv-ctd-qa/main.svg?badge_token=UugoxsuCStiAFKx-xd334Q)](https://results.pre-commit.ci/latest/github/IntegralEnvision/asv-ctd-qa/main?badge_token=UugoxsuCStiAFKx-xd334Q)
[![Docker](https://github.com/IntegralEnvision/asv-ctd-qa/workflows/docker%20build/badge.svg)](https://github.com/IntegralEnvision/asv-ctd-qa/workflows/docker-build.yml)
[![GitHub Super-Linter](https://github.com/IntegralEnvision/asv-ctd-qa/workflows/lint%20code%20base/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Convert CTD data files from ASCII to NetCDF and flag observations using [QARTOD](https://ioos.noaa.gov/project/qartod/) recommended testing implemented by [ioos_qc](https://github.com/ioos/ioos_qc).

## To-Do

- [ ] Comply with [CF Conventions](http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html)
  - [ ] Convert units: [cf-units](https://pypi.org/project/cf-units/). Udunits-2 database [xml files here](https://github.com/Unidata/UDUNITS-2/tree/master/lib)
  - [ ] [Standard names (with search)](https://cfconventions.org/Data/cf-standard-names/76/build/cf-standard-name-table.html)
  - [ ] [Standard names](http://cfconventions.org/Data/cf-standard-names/current/src/cf-standard-name-table.xml)
- [ ] Validate [config.json](./config.json) - either on a per-run basis or on update.
- [ ] Run NetCDF [compliance checker](https://github.com/ioos/compliance-checker)
  - Latest compliance check: [compliance_check.txt](./compliance_check.txt)

## Setup

### Environment

- Using [Python](https://www.python.org/downloads/release/python-3100/) version 3.9+

- Setup virtual environment

  ```bash
  python -m venv .venv
  ```

- Activate environment

  ```bash
  source ./.venv/bin/activate
  ```

- Install Dependencies

  ```bash
  pip install -r requirements.txt
  ```

  > [ioos_qc](https://github.com/ioos/ioos_qc) was cloned to the current repository on version 2.0.1 and modified for project requirements. See the [notes](#notes) section for specifics.

  - As of **2022-06-20** a few extra steps are needed to use the Python NetCDF4 module on Apple Silicon.

    ```bash
    brew install hdf5 netcdf
    git clone https://github.com/Unidata/netcdf4-python.git
    HDF5_DIR=$(brew --prefix hdf5) pip install --no-cache-dir ./netcdf4-python
    rm -r ./netcdf4-python
    ```

- Install pre-commit hooks

  ```bash
  python -m pre_commit install --install-hooks
  ```

### Configure IOOS QC Test Parameters

IOOS QC test configurations are defined in [global.json](global.json). Below is an example QC configuration for the parameter `pressure`.

```json
"pressure": {              // Standard CF-safe parameter name
  "argo": {                         // Tests based on the ARGO QC manual
    "pressure_increasing_test": null    // Check if pressure does not monotonically increase
  },
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
    "attenuated_signal_test": {         // Check for near-flat-line conditions using a range or standard deviation.
      "suspect_threshold": 5,               // Any calculated value below this amount will be flagged as SUSPECT. In observations units.
      "fail_threshold": 3,                  // Any calculated values below this amount will be flagged as FAIL. In observations units.
      "test_period": 15,                    // Length of time to test over in seconds [optional]. Otherwise, will test against entire `inp`.
      "min_obs": null,                      // Minimum number of observations in window required to calculate a result [optional]. Otherwise, test will start at beginning of time series. Note: you can specify either `min_obs` or `min_period`, but not both.
      "min_period": null,                   // Minimum number of seconds in test_period required to calculate a result [optional]. Otherwise, test will start at beginning of time series. Note: you can specify either `min_obs` or `min_period`, but not both.
      "check_type": "range"                 // Either 'std' (default) or 'range', depending on the type of check you wish to perform.
    },
    "density_inversion_test": {         // With few exceptions, potential water density will increase with increasing pressure. When vertical profile data is obtained, this test is used to flag as failed T, C, and SP observations, which yield densities that do not sufficiently increase with pressure. A small operator-selected density threshold (DT) allows for micro-turbulent exceptions.
      "suspect_threshold": 3,               // A float value representing a maximum potential density(or sigma0) variation to be tolerated, downward density variation exceeding this will be flagged as SUSPECT.
      "fail_threshold": 5                   // A float value representing a maximum potential density(or sigma0) variation to be tolerated, downward density variation exceeding this will be flagged as FAIL.
    },
    "climatology_test": {               // Checks that values are within reasonable range bounds and flags as SUSPECT.
      "suspect_span": [20, 35],              // (optional) 2-tuple range of valid values. This is passed in as the fail_span to the gross_range_test.
      "fail_span": [-5, 35],                  // 2-tuple range of valid values. This is passed in as the suspect_span to the gross_range test.
      "zspan": [0, 100]
    }
  }
}
```

### Run QC Checks

  ```bash
usage: asv_ctd_qa.py [-h] [-p] [-v] config input_file header_rows output_dir

Evaulate a data file of ASV CTD readings and apply quality assurence checks following QARTOD methods and assigning data quality flags as appropriate.
Transform results into NetCDF format following IC standards.

positional arguments:
  config                Configuration JSON file.
  input_file            Path to the input sensor data file.
  output_dir            Directory for output files.

options:
  -h, --help            show this help message and exit
  -l [LOG], --log [LOG]
                        Path to a log file for script level logging
  -p, --plot            Create an HTML file containing plots of QC flags. Files are stored under a
                        subdirectory of the specified output_dir.
  -v, --verbose         Control the amount of information to display
  ```

#### Example QC Run

  ```bash
  python asv_ctd_qa.py -v config.json ./data/received/2021-09-30T15-40-11.0.txt ./data/processed
  ```

#### Example QC Run using Docker

```bash
docker build -t asv-ctd .
```

```bash
# Bind mount the codebase during development.
docker run -it --rm -v "$(pwd)":/app asv-ctd python asv_ctd_qa.py -v config.json ./data/received/2021-09-30T15-40-11.0.txt ./data/processed
```

### Plotting QC Flags

Visualize observations and QC flags together using [qc_plots.py](./qc_plots.py). [See example output here.](./data/processed/plots/2021-09-30T15-40-11.0.txt.nc/2021-09-30T15-40-11.0.txt.nc.html)

```bash
usage: qc_plots.py [-h] [-v] ncfile outdir

Create plots of ASV CTD cast data with QARTOD flags.

positional arguments:
  ncfile         Path to the input NetCDF file.
  outdir         Path to save output files.

options:
  -h, --help     show this help message and exit
  -v, --verbose  Control the amount of information to display.
```

### Dumping NetCDF Contents

Use the `netcdf` package, which can be downloaded using Homebrew. Or, use the [ncdump.py](./ncdump.py) utility.

## Notes

- `config.json` file contains acceptable ranges for each type of sensor data.
- Parameters are ingested and used as configuration settings for QARTOD checks.
- There has been discussion of developing a set of parameters to use on a seasonal basis.

### ioos_qc

- [ioos_qc](https://github.com/ioos/ioos_qc) was cloned to the current repository on version 2.0.1 and modified for project requirements. The following revisions were made:

#### qartod.py

1. [time_interval](https://github.com/IntegralEnvision/asv-ctd-qa/commit/a249dd4ee84f719696fb31ecd6eabd9edd0f6a33#diff-32c09032f00f303300ace35369debee33af51ceb355defcce878c489bdc3af6aR646) calculation. CTD collection times can be < 1 second apart causing the number of observations to appear as 0.

## Data String Configuration (tab delimited)

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

## QA/QC & Compliance Checks

- Quality control testing is implemented using the [qartod](https://ioos.github.io/ioos_qc/api/ioos_qc.html#module-ioos_qc.qartod) module of the [ioos_qc](https://github.com/ioos/ioos_qc) Python library. See also example implentations for [Gliders](https://github.com/ioos/glider-dac).
- Compliance checks can be performed on the resulting output NetCDF via a [compliance checker](https://github.com/ioos/compliance-checker)

## Help

- [QARTOD](https://ioos.noaa.gov/project/qartod/)
- [Intro to NetCDF](https://adyork.github.io/python-oceanography-lesson/17-Intro-NetCDF/index.html)
- [Writing NetCDF](https://www.earthinversion.com/utilities/Writing-NetCDF4-Data-using-Python/)
- [Dimensions](http://www.bic.mni.mcgill.ca/users/sean/Docs/netcdf/guide.txn_12.html)
- [Western Indian Ocean (WIO) Workshop](https://github.com/MathewBiddle/WIO_workshop)
