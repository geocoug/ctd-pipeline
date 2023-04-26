#!/usr/bin/env python3
# layout.py

import dash_bootstrap_components as dbc
from dash import html


def navbar() -> dbc.Navbar:
    """Create bootstrap navigation bar."""
    return dbc.Navbar(
        dbc.Container(
            [
                html.A(
                    dbc.Row(
                        [
                            dbc.Col(
                                html.Img(
                                    src="assets/integral.png",
                                    height="40px",
                                ),
                            ),
                            dbc.Col(
                                dbc.NavbarBrand(
                                    "ASV CTD Data Viewer",
                                    className="ms-2",
                                ),
                            ),
                        ],
                        align="center",
                        className="g-0",
                    ),
                    href="/",
                    style={"textDecoration": "none"},
                ),
            ],
            fluid=True,
        ),
        color="dark",
        dark=True,
        className="navbar navbar-expand-lg",
    )


def grid() -> list:
    return [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Row(
                            html.Div(
                                [
                                    dbc.Row(
                                        options(),
                                    ),
                                ],
                                className="p-3",
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
                        ),
                    ],
                    className="col-xl-4 col-lg-5 col-12",
                    style={
                        "height": "100vh",
                        "left": 0,
                        "bottom": 0,
                        "padding": "0 2rem 0 2rem",
                        "background-color": "#F5F5F5",
                    },
                ),
                dbc.Col(
                    [
                        dbc.Row(
                            id="3d-figure-container",
                            style={
                                "height": "400px",
                            },
                        ),
                        dbc.Row(
                            html.Hr(),
                            style={
                                "height": "50px",
                            },
                        ),
                        dbc.Row(
                            id="2d-figure-container",
                            style={
                                "height": "400px",
                            },
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
                        style={"display": "none"},
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
        navbar(),
        html.Div(
            grid(),
        ),
        html.Div(id="test"),
    ]
