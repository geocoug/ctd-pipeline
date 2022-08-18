import calendar
import datetime
import json

import numpy as np
import pandas as pd

# from ioos_qc.config import ContextConfig
from ioos_qc.config import QcConfig

# from ioos_qc.results import collect_results
from ioos_qc.stores import PandasStore
from ioos_qc.streams import PandasStream


def current_time_window() -> tuple:
    """Calculate the current season time window.

    Returns a tuple with the following format:
        (season, start_time, end_time)
    """
    seasons = {
        "spring": range(3, 5),
        "summer": range(6, 8),
        "fall": range(9, 11),
    }
    today = datetime.datetime.now()
    current_season = None

    for season in seasons:
        if seasons[season].start <= today.month <= seasons[season].stop:
            current_season = season

    if not current_season:
        current_season = "winter"
        if today.month == 12:
            window_start = datetime.datetime.strptime(
                f"12 {today.year}",
                "%m %Y",
            )
            eom = calendar.monthrange(today.year + 1, 2)[1]
            window_end = datetime.datetime.strptime(
                f"2 {today.year + 1} 23:59:59",
                "%m %Y %H:%M:%S",
            )
        else:
            window_start = datetime.datetime.strptime(
                f"12 {today.year - 1}",
                "%m %Y",
            )
            eom = calendar.monthrange(today.year, 2)[1]
            window_end = datetime.datetime.strptime(
                f"2 {today.year} 23:59:59",
                "%m %Y %H:%M:%S",
            )
    else:
        window_start = datetime.datetime.strptime(
            f"{seasons[current_season].start} {today.year}",
            "%m %Y",
        )
        eom = calendar.monthrange(today.year, seasons[current_season].stop)[1]
        window_end = datetime.datetime.strptime(
            f"{eom} {seasons[current_season].stop} {today.year} 23:59:59",
            "%d %m %Y %H:%M:%S",
        )

    return (current_season, window_start, window_end)


# Load global configuration file
with open("global.json", "r", encoding="utf-8") as f:
    global_conf = json.load(f)

time_window = current_time_window()
parameters = global_conf["seasonal_params"][time_window[0]]["parameters"]

# Load in parameter file
df = pd.read_excel(parameters, sheet_name="Parameters")

# Separate dataframes per sensor
df_operator = df[df["sensor"] == "OPERATOR"]
df_turbidity = df[df["sensor"] == "TURBIDITY;FTU"]
df_pressure = df[df["sensor"] == "PRESSURE;DBAR"]
df_temp = df[df["sensor"] == "TEMPERATURE;C"]
df_do = df[df["sensor"] == "DISSOLVED OXYGEN;SAT%"]
df_altitude = df[df["sensor"] == "ALTITUDE;M"]
df_conduct = df[df["sensor"] == "CONDUCTIVITY;MS/CM"]
df_fluoro = df[df["sensor"] == "FLUOROMETER (C);UG/L"]
df_ph = df[df["sensor"] == "PH;PH"]
df_sal = df[df["sensor"] == "Calc. SALINITY; PSU"]
df_depth = df[df["sensor"] == "Calc. DEPTH;M"]


# Operator:
operator_time_inc = df_operator.loc[
    df_operator["parameter"] == "tim_inc",
    "parameter_value",
].item()
operator_samp_inc = df_operator.loc[
    df_operator["parameter"] == "samp_inc",
    "parameter_value",
].item()
operator_lat_min = df_operator.loc[
    df_operator["parameter"] == "lat_min",
    "parameter_value",
].item()
operator_lat_max = df_operator.loc[
    df_operator["parameter"] == "lat_max",
    "parameter_value",
].item()
operator_long_min = df_operator.loc[
    df_operator["parameter"] == "long_min",
    "parameter_value",
].item()
operator_long_max = df_operator.loc[
    df_operator["parameter"] == "long_max",
    "parameter_value",
].item()

# Turbidity:
turbidity_sensor_min = df_turbidity.loc[
    df_turbidity["parameter"] == "sensor_min",
    "parameter_value",
].item()
turbidity_sensor_max = df_turbidity.loc[
    df_turbidity["parameter"] == "sensor_max",
    "parameter_value",
].item()
turbidity_climate_min = df_turbidity.loc[
    df_turbidity["parameter"] == "climate_min",
    "parameter_value",
].item()
turbidity_climate_max = df_turbidity.loc[
    df_turbidity["parameter"] == "climate_max",
    "parameter_value",
].item()
turbidity_spike_high = df_turbidity.loc[
    df_turbidity["parameter"] == "spike_high",
    "parameter_value",
].item()
turbidity_spike_low = df_turbidity.loc[
    df_turbidity["parameter"] == "spike_low",
    "parameter_value",
].item()
turbidity_n_deviation = df_turbidity.loc[
    df_turbidity["parameter"] == "n_deviation",
    "parameter_value",
].item()
turbidity_t_deviation = df_turbidity.loc[
    df_turbidity["parameter"] == "t_deviation",
    "parameter_value",
].item()
turbidity_dev_interval = df_turbidity.loc[
    df_turbidity["parameter"] == "dev_interval",
    "parameter_value",
].item()
turbidity_rep_cnt_fail = df_turbidity.loc[
    df_turbidity["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
turbidity_rep_cnt_susp = df_turbidity.loc[
    df_turbidity["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
turbidity_eps = df_turbidity.loc[
    df_turbidity["parameter"] == "eps",
    "parameter_value",
].item()

# Pressure:
pressure_sensor_min = df_pressure.loc[
    df_pressure["parameter"] == "sensor_min",
    "parameter_value",
].item()
pressure_sensor_max = df_pressure.loc[
    df_pressure["parameter"] == "sensor_max",
    "parameter_value",
].item()
pressure_climate_min = df_pressure.loc[
    df_pressure["parameter"] == "climate_min",
    "parameter_value",
].item()
pressure_climate_max = df_pressure.loc[
    df_pressure["parameter"] == "climate_max",
    "parameter_value",
].item()
pressure_spike_high = df_pressure.loc[
    df_pressure["parameter"] == "spike_high",
    "parameter_value",
].item()
pressure_spike_low = df_pressure.loc[
    df_pressure["parameter"] == "spike_low",
    "parameter_value",
].item()
pressure_n_deviation = df_pressure.loc[
    df_pressure["parameter"] == "n_deviation",
    "parameter_value",
].item()
pressure_t_deviation = df_pressure.loc[
    df_pressure["parameter"] == "t_deviation",
    "parameter_value",
].item()
pressure_dev_interval = df_pressure.loc[
    df_pressure["parameter"] == "dev_interval",
    "parameter_value",
].item()
pressure_rep_cnt_fail = df_pressure.loc[
    df_pressure["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
pressure_rep_cnt_susp = df_pressure.loc[
    df_pressure["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
pressure_eps = df_pressure.loc[
    df_pressure["parameter"] == "eps",
    "parameter_value",
].item()

# Temperature:
temp_sensor_min = df_temp.loc[
    df_temp["parameter"] == "sensor_min",
    "parameter_value",
].item()
temp_sensor_max = df_temp.loc[
    df_temp["parameter"] == "sensor_max",
    "parameter_value",
].item()
temp_climate_min = df_temp.loc[
    df_temp["parameter"] == "climate_min",
    "parameter_value",
].item()
temp_climate_max = df_temp.loc[
    df_temp["parameter"] == "climate_max",
    "parameter_value",
].item()
temp_spike_high = df_temp.loc[
    df_temp["parameter"] == "spike_high",
    "parameter_value",
].item()
temp_spike_low = df_temp.loc[
    df_temp["parameter"] == "spike_low",
    "parameter_value",
].item()
temp_n_deviation = df_temp.loc[
    df_temp["parameter"] == "n_deviation",
    "parameter_value",
].item()
temp_t_deviation = df_temp.loc[
    df_temp["parameter"] == "t_deviation",
    "parameter_value",
].item()
temp_dev_interval = df_temp.loc[
    df_temp["parameter"] == "dev_interval",
    "parameter_value",
].item()
temp_rep_cnt_fail = df_temp.loc[
    df_temp["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
temp_rep_cnt_susp = df_temp.loc[
    df_temp["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
temp_eps = df_temp.loc[df_temp["parameter"] == "eps", "parameter_value"].item()

# Dissolved Oxygen:
do_sensor_min = df_do.loc[df_do["parameter"] == "sensor_min", "parameter_value"].item()
do_sensor_max = df_do.loc[df_do["parameter"] == "sensor_max", "parameter_value"].item()
do_climate_min = df_do.loc[
    df_do["parameter"] == "climate_min",
    "parameter_value",
].item()
do_climate_max = df_do.loc[
    df_do["parameter"] == "climate_max",
    "parameter_value",
].item()
do_spike_high = df_do.loc[df_do["parameter"] == "spike_high", "parameter_value"].item()
do_spike_low = df_do.loc[df_do["parameter"] == "spike_low", "parameter_value"].item()
do_n_deviation = df_do.loc[
    df_do["parameter"] == "n_deviation",
    "parameter_value",
].item()
do_t_deviation = df_do.loc[
    df_do["parameter"] == "t_deviation",
    "parameter_value",
].item()
do_dev_interval = df_do.loc[
    df_do["parameter"] == "dev_interval",
    "parameter_value",
].item()
do_rep_cnt_fail = df_do.loc[
    df_do["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
do_rep_cnt_susp = df_do.loc[
    df_do["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
do_eps = df_do.loc[df_do["parameter"] == "eps", "parameter_value"].item()

# Altitude:
altitude_sensor_min = df_altitude.loc[
    df_altitude["parameter"] == "sensor_min",
    "parameter_value",
].item()
altitude_sensor_max = df_altitude.loc[
    df_altitude["parameter"] == "sensor_max",
    "parameter_value",
].item()
altitude_climate_min = df_altitude.loc[
    df_altitude["parameter"] == "climate_min",
    "parameter_value",
].item()
altitude_climate_max = df_altitude.loc[
    df_altitude["parameter"] == "climate_max",
    "parameter_value",
].item()
altitude_spike_high = df_altitude.loc[
    df_altitude["parameter"] == "spike_high",
    "parameter_value",
].item()
altitude_spike_low = df_altitude.loc[
    df_altitude["parameter"] == "spike_low",
    "parameter_value",
].item()
altitude_n_deviation = df_altitude.loc[
    df_altitude["parameter"] == "n_deviation",
    "parameter_value",
].item()
altitude_t_deviation = df_altitude.loc[
    df_altitude["parameter"] == "t_deviation",
    "parameter_value",
].item()
altitude_dev_interval = df_altitude.loc[
    df_altitude["parameter"] == "dev_interval",
    "parameter_value",
].item()
altitude_rep_cnt_fail = df_altitude.loc[
    df_altitude["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
altitude_rep_cnt_susp = df_altitude.loc[
    df_altitude["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
altitude_eps = df_altitude.loc[
    df_altitude["parameter"] == "eps",
    "parameter_value",
].item()

# Conductivity:
conduct_sensor_min = df_conduct.loc[
    df_conduct["parameter"] == "sensor_min",
    "parameter_value",
].item()
conduct_sensor_max = df_conduct.loc[
    df_conduct["parameter"] == "sensor_max",
    "parameter_value",
].item()
conduct_climate_min = df_conduct.loc[
    df_conduct["parameter"] == "climate_min",
    "parameter_value",
].item()
conduct_climate_max = df_conduct.loc[
    df_conduct["parameter"] == "climate_max",
    "parameter_value",
].item()
conduct_spike_high = df_conduct.loc[
    df_conduct["parameter"] == "spike_high",
    "parameter_value",
].item()
conduct_spike_low = df_conduct.loc[
    df_conduct["parameter"] == "spike_low",
    "parameter_value",
].item()
conduct_n_deviation = df_conduct.loc[
    df_conduct["parameter"] == "n_deviation",
    "parameter_value",
].item()
conduct_t_deviation = df_conduct.loc[
    df_conduct["parameter"] == "t_deviation",
    "parameter_value",
].item()
conduct_dev_interval = df_conduct.loc[
    df_conduct["parameter"] == "dev_interval",
    "parameter_value",
].item()
conduct_rep_cnt_fail = df_conduct.loc[
    df_conduct["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
conduct_rep_cnt_susp = df_conduct.loc[
    df_conduct["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
conduct_eps = df_conduct.loc[df_conduct["parameter"] == "eps", "parameter_value"].item()

# Fluourometer: *** misspelled on sensor parameter document***
fluoro_sensor_min = df_fluoro.loc[
    df_fluoro["parameter"] == "sensor_min",
    "parameter_value",
].item()
fluoro_sensor_max = df_fluoro.loc[
    df_fluoro["parameter"] == "sensor_max",
    "parameter_value",
].item()
fluoro_climate_min = df_fluoro.loc[
    df_fluoro["parameter"] == "climate_min",
    "parameter_value",
].item()
fluoro_climate_max = df_fluoro.loc[
    df_fluoro["parameter"] == "climate_max",
    "parameter_value",
].item()
fluoro_spike_high = df_fluoro.loc[
    df_fluoro["parameter"] == "spike_high",
    "parameter_value",
].item()
fluoro_spike_low = df_fluoro.loc[
    df_fluoro["parameter"] == "spike_low",
    "parameter_value",
].item()
fluoro_n_deviation = df_fluoro.loc[
    df_fluoro["parameter"] == "n_deviation",
    "parameter_value",
].item()
fluoro_t_deviation = df_fluoro.loc[
    df_fluoro["parameter"] == "t_deviation",
    "parameter_value",
].item()
fluoro_dev_interval = df_fluoro.loc[
    df_fluoro["parameter"] == "dev_interval",
    "parameter_value",
].item()
fluoro_rep_cnt_fail = df_fluoro.loc[
    df_fluoro["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
fluoro_rep_cnt_susp = df_fluoro.loc[
    df_fluoro["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
fluoro_eps = df_fluoro.loc[df_fluoro["parameter"] == "eps", "parameter_value"].item()

# pH:
ph_sensor_min = df_ph.loc[df_ph["parameter"] == "sensor_min", "parameter_value"].item()
ph_sensor_max = df_ph.loc[df_ph["parameter"] == "sensor_max", "parameter_value"].item()
ph_climate_min = df_ph.loc[
    df_ph["parameter"] == "climate_min",
    "parameter_value",
].item()
ph_climate_max = df_ph.loc[
    df_ph["parameter"] == "climate_max",
    "parameter_value",
].item()
ph_spike_high = df_ph.loc[df_ph["parameter"] == "spike_high", "parameter_value"].item()
ph_spike_low = df_ph.loc[df_ph["parameter"] == "spike_low", "parameter_value"].item()
ph_n_deviation = df_ph.loc[
    df_ph["parameter"] == "n_deviation",
    "parameter_value",
].item()
ph_t_deviation = df_ph.loc[
    df_ph["parameter"] == "t_deviation",
    "parameter_value",
].item()
ph_dev_interval = df_ph.loc[
    df_ph["parameter"] == "dev_interval",
    "parameter_value",
].item()
ph_rep_cnt_fail = df_ph.loc[
    df_ph["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
ph_rep_cnt_susp = df_ph.loc[
    df_ph["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
ph_eps = df_ph.loc[df_ph["parameter"] == "eps", "parameter_value"].item()

# Calc. Salinity
sal_sensor_min = df_sal.loc[
    df_sal["parameter"] == "sensor_min",
    "parameter_value",
].item()
sal_sensor_max = df_sal.loc[
    df_sal["parameter"] == "sensor_max",
    "parameter_value",
].item()
sal_climate_min = df_sal.loc[
    df_sal["parameter"] == "climate_min",
    "parameter_value",
].item()
sal_climate_max = df_sal.loc[
    df_sal["parameter"] == "climate_max",
    "parameter_value",
].item()
sal_spike_high = df_sal.loc[
    df_sal["parameter"] == "spike_high",
    "parameter_value",
].item()
sal_spike_low = df_sal.loc[df_sal["parameter"] == "spike_low", "parameter_value"].item()
sal_n_deviation = df_sal.loc[
    df_sal["parameter"] == "n_deviation",
    "parameter_value",
].item()
sal_t_deviation = df_sal.loc[
    df_sal["parameter"] == "t_deviation",
    "parameter_value",
].item()
sal_dev_interval = df_sal.loc[
    df_sal["parameter"] == "dev_interval",
    "parameter_value",
].item()
sal_rep_cnt_fail = df_sal.loc[
    df_sal["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
sal_rep_cnt_susp = df_sal.loc[
    df_sal["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
sal_eps = df_sal.loc[df_sal["parameter"] == "eps", "parameter_value"].item()

# Calc. Depth
depth_sensor_min = df_depth.loc[
    df_depth["parameter"] == "sensor_min",
    "parameter_value",
].item()
depth_sensor_max = df_depth.loc[
    df_depth["parameter"] == "sensor_max",
    "parameter_value",
].item()
depth_climate_min = df_depth.loc[
    df_depth["parameter"] == "climate_min",
    "parameter_value",
].item()
depth_climate_max = df_depth.loc[
    df_depth["parameter"] == "climate_max",
    "parameter_value",
].item()
depth_spike_high = df_depth.loc[
    df_depth["parameter"] == "spike_high",
    "parameter_value",
].item()
depth_spike_low = df_depth.loc[
    df_depth["parameter"] == "spike_low",
    "parameter_value",
].item()
depth_n_deviation = df_depth.loc[
    df_depth["parameter"] == "n_deviation",
    "parameter_value",
].item()
depth_t_deviation = df_depth.loc[
    df_depth["parameter"] == "t_deviation",
    "parameter_value",
].item()
depth_dev_interval = df_depth.loc[
    df_depth["parameter"] == "dev_interval",
    "parameter_value",
].item()
depth_rep_cnt_fail = df_depth.loc[
    df_depth["parameter"] == "rep_cnt_fail",
    "parameter_value",
].item()
depth_rep_cnt_susp = df_depth.loc[
    df_depth["parameter"] == "rep_cnt_susp",
    "parameter_value",
].item()
depth_eps = df_depth.loc[df_depth["parameter"] == "eps", "parameter_value"].item()

bbox = tuple([operator_long_min, operator_lat_min, operator_long_max, operator_lat_max])

config = f"""
window:
    starting: {datetime.datetime.isoformat(time_window[1])}
    ending: {datetime.datetime.isoformat(time_window[2])}
streams:
    sea_water_turbidity:
        qartod:
            aggregate:
            gap_test:
            syntax_test:
            location_test:
                bbox: {bbox}
            gross_range_test:
                fail_span: {tuple([turbidity_sensor_min, turbidity_sensor_max])}
                suspect_span: {tuple([turbidity_climate_min, turbidity_climate_max])}
            climatology_test:
                config:
                tinp:
                vinp:
                zinp:
            spike_test:
                suspect_threshold: {turbidity_spike_low}
                fail_threshold: {turbidity_spike_high}
            rate_of_change_test:
            flat_line_test:
                suspect_threshold: {turbidity_rep_cnt_susp}
                fail_threshold: {turbidity_rep_cnt_fail}
                tolerance: {turbidity_eps}
"""


rows = 50
data_inputs = {
    "time": pd.date_range(start="01/01/2020", periods=rows, freq="D"),
    "z": 2.0,
    "lat": 36.1,
    "lon": -76.5,
    "sea_water_turbidity": np.arange(0, rows),
    "sea_water_temperature": np.arange(0, rows),
}
df = pd.DataFrame(data_inputs)

# Setup the stream
ps = PandasStream(df)

config = {
    "sea_water_turbidity": {
        "qartod": {
            "location_test": {
                "bbox": bbox,
            },
            "gross_range_test": {
                "fail_span": [turbidity_sensor_min, turbidity_sensor_max],
                "suspect_span": [turbidity_climate_min, turbidity_climate_max],
            },
        },
    },
}

c = QcConfig(config)

# Pass the run method the config to use
results = ps.run(c)
store = PandasStore(results)
store.compute_aggregate(name="turbidity_rollup_qc")
results_store = store.save(write_data=False, write_axes=False)
results_store = pd.concat([df, results_store], axis=1)

config = {
    "sea_water_temperature": {
        "qartod": {
            "location_test": {
                "bbox": bbox,
            },
            "gross_range_test": {
                "fail_span": [turbidity_sensor_min, turbidity_sensor_max],
                "suspect_span": [turbidity_climate_min, turbidity_climate_max],
            },
        },
    },
}

c = QcConfig(config)

# Pass the run method the config to use
results = ps.run(c)
store = PandasStore(results)
store.compute_aggregate(name="temperature_rollup_qc")
results_store = store.save(write_data=False, write_axes=False)
results_store = pd.concat([df, results_store], axis=1)


print(results_store)


# config = f"""
# region: null

# window:
#     starting: {datetime.datetime.isoformat(time_window[1])}
#     ending: {datetime.datetime.isoformat(time_window[2])}
# streams:
#     sea_water_turbidity:
#         qartod:
#             aggregate:
#             # Required Tests
#             # gap_test:
#             # syntax_test:
#             location_test:
#                 bbox: {bbox}
#             gross_range_test:
#                 fail_span: [{turbidity_sensor_min}, {turbidity_sensor_max}]
#                 suspect_span: [{turbidity_climate_min}, {turbidity_climate_max}]
#             # Decreasing Radiance, Irradiance, and PAR Test:
#             # Strongly Recommended Tests:
#             # Photic Zone Limit for Radiance, Irradiance, and PAR Test:
#             climatology_test:
#                 config:
#                 tinp:
#                 vinp:
#                 zinp:
#             spike_test:
#                 suspect_threshold: {turbidity_spike_low}
#                 fail_threshold: {turbidity_spike_high}
#             # rate_of_change_test:
#             #     threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, {turbidity_n_deviation})
#             flat_line_test:
#                 suspect_threshold: {turbidity_rep_cnt_susp}
#                 fail_threshold: {turbidity_rep_cnt_fail}
#                 tolerance: {turbidity_eps}
# """


# sea_water_pressure:
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [pressure_sensor_min, pressure_sensor_max]
#             suspect_span: [pressure_climate_min, pressure_climate_max]
#         # Strongly Recommended Tests
#         climatology_test:
#             config:
#             tinp:
#             vinp:
#             zinp:
#         spike_test:
#             suspect_threshold: pressure_spike_low
#             fail_threshold: pressure_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: pressure_rep_cnt_susp
#             fail_threshold: pressure_rep_cnt_fail
#             tolerance: pressure_eps
# Temperature: # "TEMPERATURE;C": "sea_water_temperature",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [temp_sensor_min, temp_sensor_max]
#             suspect_span: [temp_climate_min, temp_climate_max]
#         climatology_test:
#             vspan: [10, 13]
#             tspan: [1, 3]
#             period: month # what is the period?
#             zinp: [0, 100]
#         # Strongly Recommended Tests
#         spike_test:
#             suspect_threshold: temp_spike_low
#             fail_threshold: temp_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: temp_rep_cnt_susp
#             fail_threshold: temp_rep_cnt_fail
#             tolerance: temp_eps
# Dissolved Oxygen: # "DISSOLVED OXYGEN;SAT%": "sea_water_dissolved_oxygen",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [DO_sensor_min, DO_sensor_max]
#             suspect_span: [DO_climate_min, DO_climate_max]
#         climatology_test:
#             vspan: [10, 13]
#             tspan: [1, 3]
#             period: month
#             zinp: [0, 100]
#         # Strongly Recommended Tests
#         spike_test:
#             suspect_threshold: DO_spike_low
#             fail_threshold: DO_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: DO_rep_cnt_susp
#             fail_threshold: DO_rep_cnt_fail
#             tolerance: DO_eps
# Altitude: #"ALTITUDE;M": "altitude",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [altitude_sensor_min, altitude_sensor_max]
#             suspect_span: [altitude_climate_min, altitude_climate_max]
#         climatology_test:
#             config: 4
#             tinp: 6
#             vinp: 0.05
#             zinp:
#         # Strongly Recommended Tests
#         spike_test:
#             suspect_threshold: altitude_spike_low
#             fail_threshold: altitude_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: altitude_rep_cnt_susp
#             fail_threshold: altitude_rep_cnt_fail
#             tolerance: altitude_eps
# Conductivity: # "CONDUCTIVITY;MS/CM": "sea_water_electrical_conductivity",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox:[operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [conduct_sensor_min, conduct_sensor_max]
#             suspect_span: [conduct_climate_min, conduct_climate_max]
#         climatology_test:
#             config: 4
#             tinp: 6
#             vinp: 0.05
#             zinp:
#         # Strongly Recommended Tests
#         spike_test:
#             suspect_threshold: conduct_spike_low
#             fail_threshold: conduct_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: conduct_rep_cnt_susp
#             fail_threshold: conduct_rep_cnt_fail
#             tolerance: conduct_eps
# Fluorometer:  # "FLUOROMETER (C);UG/L": "sea_water_fluorescence",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [fluoro_sensor_min, fluoro_sensor_max]
#             suspect_span: [fluoro_climate_min, fluoro_climate_max]
#         #Decreasing Radiance, Irradiance, and PAR Test:

#         #Strongly Recommended Tests

#         #Photic Zone Limit for Radiance, Irradiance, and PAR Test:

#         climatology_test:
#             config: 4
#             tinp: 6
#             vinp: 0.05
#             zinp:
#         spike_test:
#             suspect_threshold: fluoro_spike_low
#             fail_threshold: fluoro_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: fluoro_rep_cnt_susp
#             fail_threshold: fluoro_rep_cnt_fail
#             tolerance: fluoro_eps
# pH: # "PH;PH": "sea_water_ph",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [pH_sensor_min, pH_sensor_max]
#             suspect_span: [pH_climate_min, pH_climate_max]
#         # Strongly Recommended Tests
#         climatology_test:
#             config: 4
#             tinp: 6
#             vinp: 0.05
#             zinp:
#         spike_test:
#             suspect_threshold: pH_spike_low
#             fail_threshold: pH_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: pH_rep_cnt_susp
#             fail_threshold: pH_rep_cnt_fail
#             tolerance: pH_eps
# Calc. Salinity: # "Calc. SALINITY; PSU": "sea_water_practical_salinity",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [sal_sensor_min, sal_sensor_max]
#             suspect_span: [sal_climate_min, sal_climate_max]
#         climatology_test:
#             config:
#             tinp:
#             vinp:
#             zinp:
#         # Strongly Recommended Tests
#         spike_test:
#             suspect_threshold: sal_spike_low
#             fail_threshold: sal_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: sal_rep_cnt_susp
#             fail_threshold: sal_rep_cnt_fail
#             tolerance: sal_eps
# Calc. Depth: # "Calc. DEPTH;M": "depth",
#     qartod:
#         aggregate:
#         # Required Tests
#         # gap_test:
#         # syntax_test:
#         location_test:
#             bbox: [operator_long_min, operator_lat_min, operator_long_max, operator_lat_max]
#         gross_range_test:
#             fail_span: [depth_sensor_min, depth_sensor_max]
#             suspect_span: [depth_climate_min, depth_climate_max]
#         climatology_test:
#             config: 4
#             tinp: 6
#             vinp: 0.05
#             zinp:
#         # Strongly Recommended Tests
#         spike_test:
#             suspect_threshold: depth_spike_low
#             fail_threshold: depth_spike_high
#         rate_of_change_test:
#             threshold: get_rate_of_change_threshold( ma.getdata(values[~mask]), ma.getdata(times[~mask]), time_units, n_dev,)
#         flat_line_test:
#             suspect_threshold: depth_rep_cnt_susp
#             fail_threshold: depth_rep_cnt_fail
#             tolerance: depth_eps
# """
