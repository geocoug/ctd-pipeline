#!/usr/bin/env python3
# callbacks.py

import configparser
import os
from datetime import datetime

import cf_xarray  # noqa
import dash
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import xarray as xr
from dash import dcc
from dash.dependencies import Input, Output, State

import viewer.layout as layout

config = configparser.ConfigParser()
config.read(os.path.join(os.path.dirname(__file__), "config.ini"))

if "data" in config.sections():
    if "ncfile_dir" in config["data"]:
        ncfile_dir = config["data"]["ncfile_dir"]
    else:
        raise configparser.NoOptionError("ncfile_dir", "data")
else:
    raise configparser.NoSectionError("data")


@dash.callback(
    Output("page-content", "children"),
    Input("url", "pathname"),
)
def display_page(pathname: str):  # noqa
    """Page content based on current path."""
    if pathname == "/":
        return layout.index()
    return dash.no_update


@dash.callback(
    Output("datafile-dropdown", "options"),
    Input("url", "pathname"),
)
def datafile_options(pathname: str):  # noqa
    """Page content based on current path."""
    if pathname == "/":
        return [
            {"label": file, "value": os.path.join(ncfile_dir, file)}
            for file in sorted(os.listdir(ncfile_dir))
            if (
                os.path.isfile(os.path.join(ncfile_dir, file))
                and os.path.splitext(file)[-1] == ".nc"
            )
        ]
    return dash.no_update


@dash.callback(
    Output("parameter-dropdown", "options"),
    Output("datafile-metadata", "children"),
    Input("datafile-dropdown", "value"),
)
def parameter_options(file: str):  # noqa
    if file:
        data = xr.open_dataset(file, decode_times=False)
        attributes = dbc.Table.from_dataframe(
            pd.DataFrame(
                {
                    "Attribute": [
                        attr
                        for attr in data.attrs
                        if attr not in ["season", "window_start", "window_end"]
                    ],
                    "Description": [
                        data.attrs[attr]
                        for attr in data.attrs
                        if attr not in ["season", "window_start", "window_end"]
                    ],
                },
            ),
            striped=True,
            bordered=True,
            hover=True,
            responsive=True,
            color="secondary",
        )
        return (
            sorted(
                [
                    name
                    for name in data.cf.standard_names
                    if (
                        "status_flag" not in name
                        and name not in ["depth", "time", "latitude", "longitude"]
                    )
                ],
            ),
            attributes,
        )

    return dash.no_update, None


@dash.callback(
    Output("qc-dropdown", "options"),
    Output("qc-dropdown-container", "style"),
    Input("qc-flag-switch", "value"),
    Input("parameter-dropdown", "value"),
    State("datafile-dropdown", "value"),
)
def qc_options(qc_flags: str, parameter: str, file: str):  # noqa
    if qc_flags:
        if len(qc_flags) > 0:
            if parameter:
                data = xr.open_dataset(file, decode_times=False)
                return (
                    sorted(
                        [
                            qc.replace(f"{parameter}_", "")
                            for qc in data.cf.standard_names[f"{parameter} status_flag"]
                        ],
                    ),
                    None,
                )
            return None, None
    return None, {"display": "none"}


@dash.callback(
    Output("2d-figure-container", "children"),
    Output("3d-figure-container", "children"),
    Input("datafile-dropdown", "value"),
    Input("parameter-dropdown", "value"),
    Input("qc-flag-switch", "value"),
    Input("qc-dropdown", "value"),
)
def generate_plots(datafile, parameter, include_qc, qc_type):  # noqa
    if not parameter:
        return [None], [None]

    data = xr.open_dataset(datafile, decode_times=False)
    time = [datetime.fromtimestamp(x) for x in data["time"].values]
    units = data[parameter].units
    if str(units) == "1":
        units = "#"
    obs = data[parameter].values

    # Create 3d plot
    fig_3d = px.scatter_3d(
        x=data["longitude"],
        y=data["latitude"],
        z=-data["depth"],
        color=data[parameter],
    )
    fig_3d.update_layout(
        coloraxis_colorbar=dict(
            title=units,
            # https://plotly.com/python/colorscales/
            # thicknessmode="pixels",
            # thickness=50,
            lenmode="pixels",
            len=300,
            # yanchor="top",
            # y=1,
            # ticks="outside",
            # ticksuffix=" bills",
            # dtick=5,
        ),
    )

    # Create 2d plot
    fig_2d = go.Figure()
    fig_2d.add_trace(
        go.Scatter(
            x=time,
            y=obs,
            name="obs",
            mode="markers",
            opacity=1,
            marker=dict(size=4, color="#637c8a"),
        ),
    )
    fig_2d.update_layout(
        yaxis_title=f"<b>{units}</b>",
        xaxis_title="<b>time</b>",
        margin=dict(l=10, r=10, b=10, t=10),
        legend=dict(
            # yanchor="top",
            y=0.6,
            # xanchor="right",
        ),
    )

    if include_qc:
        if len(include_qc) > 0 and qc_type:
            qc_test = data[f"{parameter}_{qc_type}"].values
            qc_pass = np.ma.masked_where(qc_test != 1, obs).filled(np.nan)
            qc_notrun = np.ma.masked_where(qc_test != 2, obs).filled(np.nan)
            qc_suspect = np.ma.masked_where(qc_test != 3, obs).filled(np.nan)
            qc_fail = np.ma.masked_where(qc_test != 4, obs).filled(np.nan)

            fig_2d.add_trace(
                go.Scatter(
                    x=time,
                    y=qc_notrun,
                    name="qc not run",
                    mode="markers",
                    opacity=0.3,
                    marker=dict(size=4, color="gray"),
                ),
            )
            fig_2d.add_trace(
                go.Scatter(
                    x=time,
                    y=qc_pass,
                    name="qc pass",
                    mode="markers",
                    opacity=0.5,
                    marker=dict(size=5, color="green"),
                ),
            )
            fig_2d.add_trace(
                go.Scatter(
                    x=time,
                    y=qc_suspect,
                    name="qc suspect",
                    mode="markers",
                    opacity=0.7,
                    marker=dict(size=5, color="orange"),
                ),
            )
            fig_2d.add_trace(
                go.Scatter(
                    x=time,
                    y=qc_fail,
                    name="qc fail",
                    mode="markers",
                    opacity=1,
                    marker=dict(size=5, color="red"),
                ),
            )

            fig_3d = go.Figure()
            fig_3d.add_trace(
                go.Scatter3d(
                    x=data["longitude"],
                    y=data["latitude"],
                    z=-data["depth"],
                    name="obs",
                    mode="markers",
                    marker=dict(size=3, color="#637c8a"),
                    opacity=1,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 2, data["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 2, data["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 2, -data["depth"]).filled(np.nan),
                    name="qc not run",
                    mode="markers",
                    marker=dict(size=6, color="gray"),
                    opacity=0.3,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 1, data["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 1, data["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 1, -data["depth"]).filled(np.nan),
                    name="qc pass",
                    mode="markers",
                    marker=dict(size=4, color="green"),
                    opacity=0.5,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 3, data["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 3, data["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 3, -data["depth"]).filled(np.nan),
                    name="qc suspect",
                    mode="markers",
                    marker=dict(size=5, color="orange"),
                    opacity=0.7,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 4, data["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 4, data["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 4, -data["depth"]).filled(np.nan),
                    name="qc fail",
                    mode="markers",
                    marker=dict(size=6, color="red"),
                    opacity=1,
                ),
            )
    fig_3d.update_layout(
        scene=dict(
            yaxis_title="<b>latitude</b>",
            xaxis_title="<b>longitude</b>",
            zaxis_title="<b>depth</b>",
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=0.5),
        ),
        legend=dict(
            yanchor="top",
            y=0.6,
            xanchor="right",
        ),
        margin=dict(l=10, r=10, b=10, t=10),
    )
    return dcc.Graph(figure=fig_2d, id="2d-figure"), dcc.Graph(
        figure=fig_3d,
        id="3d-figure",
    )