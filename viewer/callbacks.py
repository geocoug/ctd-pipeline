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
from dotenv import load_dotenv

import viewer.layout as layout

APP_DIR = os.getenv("APP_DIR")
if not APP_DIR:
    APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

load_dotenv(os.path.join(APP_DIR, ".env"))

MAPBOX_TOKEN = os.getenv("MAPBOX_TOKEN")
px.set_mapbox_access_token(MAPBOX_TOKEN)
COLOR_SCALE = px.colors.sequential.Bluered


config = configparser.ConfigParser()
config.read(os.path.join(APP_DIR, "viewer/config.ini"))

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


# Set datafile dropdown options and default value
@dash.callback(
    Output("datafile-dropdown", "options"),
    Output("datafile-dropdown", "value"),
    Input("url", "pathname"),
)
def datafile_options(pathname: str):  # noqa
    """Page content based on current path."""
    if pathname == "/":
        options = [
            {"label": file, "value": os.path.join(ncfile_dir, file)}
            for file in sorted(os.listdir(ncfile_dir))
            if (
                os.path.isfile(os.path.join(ncfile_dir, file))
                and os.path.splitext(file)[-1] == ".nc"
            )
        ]
        return options, options[-1]["value"]
    return dash.no_update, dash.no_update


# Set parameter dropdown options and default value
@dash.callback(
    Output("parameter-dropdown", "options"),
    Output("parameter-dropdown", "value"),
    Output("map-container", "children"),
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

        # https://plotly.com/python/scattermapbox/
        lat = float(sum(data["latitude"].values) / len(data["latitude"].values))
        lon = float(sum(data["longitude"]) / len(data["longitude"]))
        mapfig = go.Figure(
            go.Scattermapbox(
                lat=[lat],
                lon=[lon],
                marker=go.scattermapbox.Marker(
                    size=12,
                ),
            ),
        )
        mapfig.update_layout(
            autosize=True,
            hovermode="closest",
            showlegend=False,
            mapbox=dict(
                style="light",
                accesstoken=MAPBOX_TOKEN,
                zoom=8,
                center=dict(
                    lat=lat,
                    lon=lon,
                ),
            ),
            margin=dict(l=10, r=10, b=10, t=10),
        )
        parameters = sorted(
            [
                name
                for name in data.cf.standard_names
                if (
                    "status_flag" not in name
                    and name not in ["depth", "time", "latitude", "longitude"]
                )
            ],
        )
        default_parameter = parameters[0]
        for parameter in parameters:
            if "temperature" in parameter:
                default_parameter = parameter
        return (
            parameters,
            default_parameter,
            dcc.Graph(
                figure=mapfig,
                id="map-figure",
                style={
                    "height": "180px",
                },
            ),
            attributes,
        )

    return dash.no_update, None, None, None


@dash.callback(
    Output("qc-dropdown", "options"),
    Output("qc-dropdown", "value"),
    Input("qc-flag-switch", "value"),
    Input("parameter-dropdown", "value"),
    State("datafile-dropdown", "value"),
)
def qc_options(qc_flags: str, parameter: str, file: str):  # noqa
    if qc_flags:
        if len(qc_flags) > 0:
            if parameter:
                data = xr.open_dataset(file, decode_times=False)
                qc_tests = sorted(
                    [
                        qc.replace(f"{parameter}_", "")
                        for qc in data.cf.standard_names[f"{parameter} status_flag"]
                    ],
                )
                return qc_tests, qc_tests[-1]
            return None, None
    return None, None


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
        color_continuous_scale=COLOR_SCALE,
    )
    fig_3d.update_layout(
        coloraxis_colorbar=dict(
            title=f"<b>{units}</b>",
            # https://plotly.com/python/colorscales/
            lenmode="pixels",
            len=300,
        ),
    )

    # Create 2d plot
    fig_2d = px.scatter(
        x=time,
        y=-data["depth"].values,
        color=obs,
        color_continuous_scale=COLOR_SCALE,
    )
    fig_2d.update_layout(
        yaxis_title="<b>depth</b>",
        xaxis_title="<b>time</b>",
        margin=dict(l=10, r=10, b=10, t=10),
        coloraxis_colorbar=dict(
            title=f"<b>{units}</b>",
            y=0.5,
            lenmode="pixels",
            len=300,
        ),
    )

    if include_qc:
        if len(include_qc) > 0 and qc_type:
            qc_test = data[f"{parameter}_{qc_type}"].values
            qc_pass = np.ma.masked_where(qc_test != 1, obs).filled(np.nan)
            qc_notrun = np.ma.masked_where(qc_test != 2, obs).filled(np.nan)
            qc_suspect = np.ma.masked_where(qc_test != 3, obs).filled(np.nan)
            qc_fail = np.ma.masked_where(qc_test != 4, obs).filled(np.nan)

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
                    y=0.6,
                ),
            )

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
