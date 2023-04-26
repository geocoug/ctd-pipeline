#!/usr/bin/env python3
# app.py


import dash
import dash_bootstrap_components as dbc
from dash import dcc, html

import viewer.callbacks as callbacks  # noqa
import viewer.layout as layout  # noqa

app = dash.Dash(
    __name__,
    title="ASV CTD Data Viewer",
    assets_folder="assets",
    external_stylesheets=[
        dbc.themes.BOOTSTRAP,
        dbc.icons.FONT_AWESOME,
    ],
    suppress_callback_exceptions=True,  # Disable when creating objects using callbacks
    meta_tags=[
        {
            "name": "viewport",
            "content": "width=device-width, initial-scale=1",
        },
    ],
)

app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        html.Div(id="page-content"),
    ],
    id="container",
)

dashboard = app.server


if __name__ == "__main__":
    app.run(debug=True)
