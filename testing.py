import calendar
import datetime
import json
from numbers import Real

import netCDF4 as nc
import numpy as np
import pandas as pd

import qc_plots
from ioos_qc import qartod, utils
from ioos_qc.config import QcConfig
from ioos_qc.stores import PandasStore
from ioos_qc.streams import PandasStream

# from qartod_config import global_config, qartod_config, standard_names

N = Real


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


with open("global.json", "r", encoding="utf-8") as f:
    config = json.load(f)

infile = "./data/received/2021-09-30T15-40-11.0.txt"
df = pd.read_csv(infile, sep="\t")
df.rename(columns=config["standard_names"], inplace=True)

df["time"] = df["time"].astype("datetime64")
df["time"] = df["time"] - datetime.datetime(1970, 1, 1, 0, 0)
df["time"] = df["time"].dt.total_seconds()

for parameter in config["qc_tests"]["parameters"]:
    parameter_config = config["qc_tests"]["parameters"][parameter]
    if "valid_range_test" in parameter_config["axds"]:
        parameter_config["axds"]["valid_range_test"]["valid_span"] = tuple(
            np.array(parameter_config["axds"]["valid_range_test"]["valid_span"]),
        )
    if "density_inversion_test" in parameter_config["qartod"]:
        parameter_config["qartod"]["density_inversion_test"]["zinp"] = df[
            "depth"
        ].to_numpy()
    if "climatology_test" in parameter_config["qartod"]:
        parameter_config["qartod"]["climatology_test"]["zinp"] = df["depth"].to_numpy()
        season, window_start, window_end = current_time_window()
        cc = qartod.ClimatologyConfig()
        cc.add(
            tspan=tuple([window_start, window_end]),
            fspan=tuple(
                np.array(parameter_config["qartod"]["climatology_test"]["fail_span"]),
            ),
            vspan=tuple(
                np.array(
                    parameter_config["qartod"]["climatology_test"]["suspect_span"],
                ),
            ),
        )
        parameter_config["qartod"]["climatology_test"]["config"] = cc
    # Add global tests
    parameter_config["qartod"].update(config["qc_tests"]["global"]["qartod"])


ps = PandasStream(df)
qc = QcConfig(config["qc_tests"]["parameters"])
qc_results = ps.run(qc)
store = PandasStore(qc_results)
store.compute_aggregate()
results_store = store.save(write_data=False, write_axes=False)
results = pd.concat([df, results_store], axis=1)
# print(results.head())
# print(results.columns)


# Don't need to include this as a column, but should probably go into the global/metadata
# whether it passed or not (returns True | False)
x = utils.check_timestamps(
    df["time"].to_numpy(),
    **config["qc_tests"]["global"]["utils"]["check_timestamps"],
)
# print(x)


rootgrp = nc.Dataset("temp.nc", "w", format="NETCDF4")
rootgrp.createDimension("time", None)
time = rootgrp.createVariable("time", np.dtype("float64").char, ("time",))
time.standard_name = "time"
time.units = "seconds since 1970-01-01 00:00:00"
time[:] = results["time"]

lat = rootgrp.createVariable("lat", np.dtype("float64").char, ("time",))
lat.standard_name = "latitude"
lat.units = "degrees_north"
lat[:] = results["lat"]

lon = rootgrp.createVariable("lon", np.dtype("float64").char, ("time",))
lon.standard_name = "longitude"
lon.units = "degrees_north"
lon[:] = results["lon"]

exclude_cols = ["time", "lat", "lon", "qartod_rollup"]
for col in results.columns:
    parent = col
    if col not in exclude_cols:
        if "test" in col:
            if "qartod" in col:
                parent = col.split("_qartod")[0]
            elif "argo" in col:
                parent = col.split("_argo")[0]
            elif "axds" in col:
                parent = col.split("_axds")[0]
        var = rootgrp.createVariable(col, np.dtype("float64").char, ("time",))
        var.standard_name = f"{parent} status_flag"
        var.long_name = col
        var.units = "idk"
        var[:] = results[col]


qc_plots.generate_plots("temp.nc", "test", verbose=False)
