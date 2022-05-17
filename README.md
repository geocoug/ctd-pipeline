# asv-ctd-qa

[![GitHub Super-Linter](https://github.com/IntegralEnvision/asv-ctd-qa/workflows/lint%20code%20base/badge.svg)](https://github.com/marketplace/actions/super-linter)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Convert CTD data files from ASCII to NetCDF and flag observations using [QARTOD](https://ioos.noaa.gov/project/qartod/) recommended testing implemented by [ioos_qc](https://github.com/ioos/ioos_qc).

## To-Do

1. Execute QARTOD checks using the [ioos_qc](https://github.com/ioos/ioos_qc) library
   - See example moduel usage [here](https://github.com/ioos/glider-dac)
1. Validate [Sensor Parameters](./parameters/Sensor_Parameters.xlsx) - either on a per-run basis or on update.
1. Collaborate with GCOOS on delivery format
1. Validate output - run manual test against generated output
1. Run NetCDF [compliance checker](https://github.com/ioos/compliance-checker)
1. Integrate database archiving.

## Usage

```bash
python main.py -v ./data/received/2021-09-30T15-40-11.0.txt 0 ./parameters/Sensor_Parameters.xlsx ./data/processed logs
```

## Notes

- Operator defined [Sensor Parameters](./parameters/Sensor_Parameters.xlsx) are defined for each cruise.
- Parameters file contains acceptable ranges for the various data types.
- Parameters are ingested and used as configuration settings for QARTOD checks.
- There has been discussion of developing a set of parameters to use on a seasonal basis.

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
