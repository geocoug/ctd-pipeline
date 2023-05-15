#!/usr/bin/env python3
# callbacks.py

import configparser

# import json
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


# def ncfile_metadata() -> dict:
#     metadata = {}
#     for file in os.listdir(ncfile_dir):
#         fpath = os.path.join(ncfile_dir, file)
#         if os.path.splitext(fpath)[-1] != ".nc":
#             continue
#         data = xr.open_dataset(os.path.join(ncfile_dir, file), decode_times=False)
#         file_metadata = {
#             attr: data.attrs[attr]
#             for attr in data.attrs
#             if attr not in ["season", "window_start", "window_end"]
#         }
#         file_metadata.update(
#             {
#                 "Start Time": datetime.fromtimestamp(data["time"].values[0]).isoformat(  # noqa
#                     sep=" ",
#                 ),
#                 "End Time": datetime.fromtimestamp(data["time"].values[-1]).isoformat(
#                     sep=" ",
#                 ),
#                 "Elapsed Time": f"""{(datetime.fromtimestamp(data['time'].values[-1])
#             - datetime.fromtimestamp(data['time'].values[0])).seconds} seconds""",
#                 "Max Depth": f"{max(data['depth'].values)} [{data['depth'].units}]",
#             },
#         )
#         metadata.update(
#             {
#                 file: {
#                     "metadata": {
#                         "Attribute": list(file_metadata.keys()),
#                         "Description": [file_metadata[key] for key in file_metadata],
#                     },
#                     "lat": float(
#                         sum(data["latitude"].values) / len(data["latitude"].values),
#                     ),
#                     "lon": float(
#                         sum(data["longitude"].values) / len(data["longitude"].values),
#                     ),
#                     **file_metadata,
#                 },
#             },
#         )
#     return metadata


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
    # Output("metadata", "data"),
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
        return options, options[-1]["value"]  # , json.dumps(ncfile_metadata())
    return dash.no_update, dash.no_update  # , dash.no_update


# Set parameter dropdown options and default value
@dash.callback(
    Output("parameter-dropdown", "options"),
    Output("parameter-dropdown", "value"),
    Output("map-container", "children"),
    Output("datafile-metadata", "children"),
    Output("dataset", "data"),
    Input("datafile-dropdown", "value"),
    # State("metadata", "data"),
)
def parameter_options(file: str):  # , metadata):  # noqa
    if file:
        # metadata = json.loads(metadata)
        sel_file = os.path.basename(file)
        data = xr.open_dataset(file, decode_times=False)

        metadata = {
            attr: data.attrs[attr]
            for attr in data.attrs
            if attr not in ["season", "window_start", "window_end"]
        }
        cast_attrs = {
            "Start Time": datetime.fromtimestamp(data["time"].values[0]).isoformat(
                sep=" ",
            ),
            "End Time": datetime.fromtimestamp(data["time"].values[-1]).isoformat(
                sep=" ",
            ),
            "Elapsed Time": f"""{(datetime.fromtimestamp(data['time'].values[-1])
            - datetime.fromtimestamp(data['time'].values[0])).seconds} seconds""",
            "Max Depth": f"{max(data['depth'].values)} [{data['depth'].units}]",
        }
        metadata.update(cast_attrs)

        df_attrs = {
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
        }
        df_attrs["Attribute"].append(list(cast_attrs.keys()))
        df_attrs["Description"].append(list(cast_attrs.values()))

        metadata_table = dbc.Table.from_dataframe(
            # pd.DataFrame(metadata[sel_file]["metadata"]),
            pd.DataFrame(df_attrs),
            striped=True,
            bordered=True,
            hover=True,
            responsive=True,
            className="table-sm table-responsive",
            style={
                "overflowX": "scroll",
                "overflowY": "scroll",
            },
        )
        # https://plotly.com/python/scattermapbox/
        lat = float(sum(data["latitude"].values) / len(data["latitude"].values))
        lon = float(sum(data["longitude"].values) / len(data["longitude"].values))
        mapfig = go.Figure(
            go.Scattermapbox(
                lat=[lat],
                lon=[lon],
                marker=go.scattermapbox.Marker(
                    size=12,
                ),
                hoverinfo=None,
                hovertext=f"""<b>File</b>: {sel_file}<br><b>Start Time</b>: {metadata["Start Time"]}<br><b>End Time</b>: {metadata["End Time"]}<br><b>Elapsed Time</b>: {metadata["Elapsed Time"]}<br><b>Max Depth</b>: {metadata["Max Depth"]}""",  # noqa
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
                    # lat=metadata[sel_file]["lat"],
                    # lon=metadata[sel_file]["lon"],
                    lat=lat,
                    lon=lon,
                ),
            ),
            margin=dict(l=5, r=5, b=5, t=5),
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
            metadata_table,
            data.to_dict(),
        )

    return dash.no_update, None, None, None, None


@dash.callback(
    Output("qc-dropdown", "options"),
    Output("qc-dropdown", "value"),
    Input("qc-flag-switch", "value"),
    Input("parameter-dropdown", "value"),
    State("dataset", "data"),
)
def qc_options(qc_flags: str, parameter: str, data: dict):  # noqa
    if qc_flags:
        if len(qc_flags) > 0:
            if parameter:
                qc_tests = sorted(
                    [
                        qc.replace(f"{parameter}_", "")
                        for qc in data["data_vars"]
                        if parameter in qc and qc != parameter
                    ],
                )
                return qc_tests, qc_tests[-1]
    return None, None


@dash.callback(
    Output("2d-figure-container", "children"),
    Output("3d-figure-container", "children"),
    Input("datafile-dropdown", "value"),
    Input("parameter-dropdown", "value"),
    Input("qc-flag-switch", "value"),
    Input("qc-dropdown", "value"),
    State("dataset", "data"),
)
def generate_plots(datafile, parameter, include_qc, qc_type, data):  # noqa
    if not parameter:
        return [None], [None]
    time = [datetime.fromtimestamp(x) for x in data["coords"]["time"]["data"]]
    depth_unit = data["data_vars"]["depth"]["attrs"]["units"]
    units = data["data_vars"][parameter]["attrs"]["units"]
    if str(units) == "1":
        units = "#"
    obs = data["data_vars"][parameter]["data"]

    df = pd.DataFrame(
        {
            "time": time,
            "longitude": data["data_vars"]["longitude"]["data"],
            "latitude": data["data_vars"]["latitude"]["data"],
            "depth": -np.array(data["data_vars"]["depth"]["data"]),
            "obs": obs,
            "color": obs,
            "units": units,
        },
    )
    # Create 3d plot
    fig_3d = px.scatter_3d(
        df,
        x="longitude",
        y="latitude",
        z="depth",
        color="color",
        color_continuous_scale=COLOR_SCALE,
        hover_data={"obs": True, "color": False, "units": True},
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
        df,
        x="time",
        y="depth",
        color="obs",
        color_continuous_scale=COLOR_SCALE,
    )
    fig_2d.update_layout(
        yaxis_title=f"<b>depth [{depth_unit}]</b>",
        xaxis_title="<b>time</b>",
        margin=dict(l=10, r=10, b=10, t=30),
        coloraxis_colorbar=dict(
            title=f"<b>{units}</b>",
            y=0.5,
            lenmode="pixels",
            len=300,
        ),
    )

    if include_qc:
        if len(include_qc) > 0 and qc_type:
            qc_test = np.array(data["data_vars"][f"{parameter}_{qc_type}"]["data"])
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
                    x=df["longitude"],
                    y=df["latitude"],
                    z=df["depth"],
                    name="obs",
                    mode="markers",
                    marker=dict(size=3, color="#637c8a"),
                    opacity=1,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 2, df["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 2, df["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 2, df["depth"]).filled(np.nan),
                    name="qc not run",
                    mode="markers",
                    marker=dict(size=6, color="gray"),
                    opacity=0.3,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 1, df["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 1, df["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 1, df["depth"]).filled(np.nan),
                    name="qc pass",
                    mode="markers",
                    marker=dict(size=4, color="green"),
                    opacity=0.5,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 3, df["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 3, df["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 3, df["depth"]).filled(np.nan),
                    name="qc suspect",
                    mode="markers",
                    marker=dict(size=5, color="orange"),
                    opacity=0.7,
                ),
            )
            fig_3d.add_trace(
                go.Scatter3d(
                    x=np.ma.masked_where(qc_test != 4, df["longitude"]).filled(
                        np.nan,
                    ),
                    y=np.ma.masked_where(qc_test != 4, df["latitude"]).filled(np.nan),
                    z=np.ma.masked_where(qc_test != 4, df["depth"]).filled(np.nan),
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
            zaxis_title=f"<b>depth [{depth_unit}]</b>",
            aspectmode="manual",
            aspectratio=dict(x=1, y=1, z=0.5),
            camera=dict(
                eye=dict(
                    x=2,
                    y=2,
                    z=0.5,
                ),
            ),
        ),
        legend=dict(
            yanchor="top",
            y=0.6,
            xanchor="right",
        ),
        margin=dict(l=10, r=10, b=10, t=10),
    )

    fig_2d.update_traces(
        hovertemplate="<b>Time</b>: %{x}<br><b>Value</b>: %{y} " + f"[{units}]",
    )
    return dcc.Graph(figure=fig_2d, id="2d-figure", clear_on_unhover=True), dcc.Graph(
        figure=fig_3d,
        id="3d-figure",
        clear_on_unhover=True,
    )


# @dash.callback(
#     Output("2d-figure", "figure"),
#     Output("3d-figure", "figure"),
#     Input("qc-flag-switch", "value"),
#     Input("2d-figure", "hoverData"),
#     Input("3d-figure", "hoverData"),
#     State("2d-figure", "figure"),
#     State("3d-figure", "figure"),
# )
# def hover_animations(include_qc, hover_2d, hover_3d, state_2d, state_3d):
#     n_traces = 1
#     if include_qc and len(include_qc) > 0:
#         n_traces = 5
#     state_2d["data"] = state_2d["data"][:n_traces]
#     state_3d["data"] = state_3d["data"][:n_traces]
#     if hover_2d:
#         idx = hover_2d["points"][0]["pointNumber"]
#         state_3d["data"].append(
#             {
#                 "x": [state_3d["data"][0]["x"][idx]],
#                 "y": [state_3d["data"][0]["y"][idx]],
#                 "z": [state_3d["data"][0]["z"][idx]],
#                 "mode": "markers",
#                 "marker": {
#                     "color": "#00fcff",
#                     "size": 12,
#                 },
#                 "showlegend": False,
#                 "type": "scatter3d",
#             },
#         )
#         return dash.no_update, state_3d
#     if hover_3d:
#         idx = hover_3d["points"][0]["pointNumber"]
#         state_2d["data"].append(
#             {
#                 "x": [state_2d["data"][0]["x"][idx]],
#                 "y": [state_2d["data"][0]["y"][idx]],
#                 "mode": "markers",
#                 "marker": {
#                     "color": "#00fcff",
#                     "size": 12,
#                 },
#                 "showlegend": False,
#                 "type": "scatter",
#             },
#         )
#         return state_2d, dash.no_update
#     return state_2d, state_3d
