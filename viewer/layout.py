#!/usr/bin/env python3
# layout.py

import dash_bootstrap_components as dbc
from dash import html


def grid() -> list:
    return [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    html.H3(
                                        "ASV CTD Data Viewer",
                                        className="p-2 mx-auto my-auto",
                                        style={"color": "#c8d3d5"},
                                    ),
                                    className="col justify-content-center text-center",
                                ),
                                dbc.Col(
                                    html.Img(
                                        src="assets/integral.png",
                                        className="mx-auto",
                                        style={
                                            "width": "75px",
                                            "height": "50px",
                                        },
                                    ),
                                    className="col-md-3 col-4 mx-auto my-auto",
                                ),
                            ],
                            style={"background-color": "#2c353a"},
                        ),
                        dbc.Row(
                            html.Div(
                                [
                                    dbc.Row(
                                        options(),
                                    ),
                                ],
                                className="p-1",
                            ),
                        ),
                        dbc.Row(
                            html.Div(
                                id="map-container",
                            ),
                        ),
                        dbc.Row(
                            html.Div(
                                [
                                    dbc.Row(
                                        html.Div(
                                            id="datafile-metadata",
                                            className="mx-auto my-2",
                                        ),
                                    ),
                                ],
                            ),
                            className="my-auto mt-2",
                        ),
                    ],
                    className="col-xl-4 col-lg-5 col-12",
                    style={
                        "height": "100vh",
                        "left": 0,
                        "bottom": 0,
                        "background-color": "#F5F5F5",
                    },
                ),
                dbc.Col(
                    [
                        dbc.Row(
                            id="3d-figure-container",
                        ),
                        dbc.Row(
                            html.Hr(),
                        ),
                        dbc.Row(
                            id="2d-figure-container",
                        ),
                    ],
                    className="col-xl-8 col-lg-7 col-12",
                ),
            ],
            id="map-row",
        ),
    ]


def options() -> list:
    return [
        html.Div(
            [
                dbc.InputGroup(
                    [
                        dbc.InputGroupText("Data File"),
                        dbc.Select(
                            options=[],
                            id="datafile-dropdown",
                        ),
                    ],
                    className="p-2",
                ),
                dbc.InputGroup(
                    [
                        dbc.InputGroupText("Parameter"),
                        dbc.Select(
                            options=[],
                            id="parameter-dropdown",
                        ),
                    ],
                    className="p-2",
                ),
                dbc.InputGroup(
                    [
                        dbc.Checklist(
                            options=[
                                {"label": "Plot QC Flags", "value": "on"},
                            ],
                            id="qc-flag-switch",
                            switch=True,
                            inline=True,
                            className="mx-auto text-center",
                        ),
                    ],
                    className="p-2",
                ),
                html.Div(
                    dbc.InputGroup(
                        [
                            dbc.InputGroupText("QC Flag"),
                            dbc.Select(
                                id="qc-dropdown",
                            ),
                        ],
                        id="qc-dropdown-container",
                        className="p-2",
                        # style={"display": "none"},
                    ),
                ),
            ],
            id="options-menu",
            className="my-2 mx-auto",
        ),
    ]


def index() -> list:
    """Create index page content."""
    return [
        dbc.Container(
            [
                html.Div(id="test"),
                html.Div(grid()),
            ],
            fluid=True,
        ),
    ]
