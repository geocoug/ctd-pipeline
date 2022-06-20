# asv-ctd-qa

[![GitHub Super-Linter](https://github.com/IntegralEnvision/asv-ctd-qa/workflows/lint%20code%20base/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Convert CTD data files from ASCII to NetCDF and flag observations using [QARTOD](https://ioos.noaa.gov/project/qartod/) recommended testing implemented by [ioos_qc](https://github.com/ioos/ioos_qc).

## To-Do

- [ ] Unit testing
- [ ] Execute QARTOD checks using the [ioos_qc](https://github.com/ioos/ioos_qc) library (see example module usage [here](https://github.com/ioos/glider-dac))
- [ ] Validate [Sensor Parameters](./parameters/Sensor_Parameters.xlsx) - either on a per-run basis or on update.
- [ ] Convert parameters to YAML?
- [ ] Collaborate with GCOOS on delivery format
- [ ] Validate output - run manual test against generated output
- [ ] Run NetCDF [compliance checker](https://github.com/ioos/compliance-checker)
- [ ] Integrate database archiving.

## Setup

- Using [Python](https://www.python.org/downloads/release/python-3100/) version 3.10+

- Setup Python virtual environment

  ```bash
  python -m venv env
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

- Run QC

  ```bash
  python main.py -v ./data/received/2021-09-30T15-40-11.0.txt 0 ./parameters/Sensor_Parameters.xlsx ./data/processed logs
  ```

## Notes

- Operator defined [Sensor Parameters](./parameters/Sensor_Parameters.xlsx) are defined for each cruise.
- Parameters file contains acceptable ranges for the various data types.
- Parameters are ingested and used as configuration settings for QARTOD checks.
- There has been discussion of developing a set of parameters to use on a seasonal basis.

### ioos_qc

- [ioos_qc](https://github.com/ioos/ioos_qc) was cloned to the current repository on version 2.0.1 and modified for project requirements. The following revisions were made:

  **qartod.py**

  - [time_interval](https://github.com/IntegralEnvision/asv-ctd-qa/commit/a249dd4ee84f719696fb31ecd6eabd9edd0f6a33#diff-32c09032f00f303300ace35369debee33af51ceb355defcce878c489bdc3af6aR646) calculation. CTD collection times can be < 1 second apart causing the number of observations to appear as 0.

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

Quality control testing is implemented using the [qartod](https://ioos.github.io/ioos_qc/api/ioos_qc.html#module-ioos_qc.qartod) module of the [ioos_qc](https://github.com/ioos/ioos_qc) Python library. See also example implentations for [Gliders](https://github.com/ioos/glider-dac).

Compliance checks can be performed on the resulting output NetCDF via a [compliance checker](https://github.com/ioos/compliance-checker)

## Help

- [QARTOD](https://ioos.noaa.gov/project/qartod/)
- [Intro to NetCDF](https://adyork.github.io/python-oceanography-lesson/17-Intro-NetCDF/index.html)
- [Writing NetCDF](https://www.earthinversion.com/utilities/Writing-NetCDF4-Data-using-Python/)
- [Dimensions](http://www.bic.mni.mcgill.ca/users/sean/Docs/netcdf/guide.txn_12.html)
