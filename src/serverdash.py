#!/usr/bin/env python3

# Author: Sean Maden
# 
# Dashboard script for recountmethylation server.
#
#

import datetime, dash, plotly, os, sys, re
import pandas as pd; import dash_core_components as dcc
import dash_html_components as html; import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output; import plotly.express as px
import plotly.graph_objects as go; from plotly.subplots import make_subplots
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from report import idats_report, soft_report

# GLOBALS
# define globals
NUM_SEC = 60
# new query
ddidat = idats_report(); ddsoft = soft_report(ddidat)
ddsoft_lkeys = list(ddsoft.keys())
# main data for line plots by time
data_all = pd.DataFrame({"num_idat" : [0], "num_gsesoft" : [0],
    "num_gsmsoft" : [0], "time" : [0]})
# file counts data
dfcount = pd.DataFrame({"file_type": [
        "IDAT", "GSE SOFT", "GSM SOFT", "GSM SOFT & IDAT",
        "IDAT", "IDAT", "IDAT", "GSE SOFT", "GSM SOFT", "GSM JSON"
    ],
    "count_type" : [
        "Unique GSM ID", "Unique GSE ID", "Unique GSM ID", "Unique GSM ID",
        "Ext. 'idat.gz'", "Ext. `idat`", "Type hlink", "Ext. `soft.gz`", 
        "Ext. `soft`", "Ext. `json.filt`"
    ],
    "counts" : [
        ddidat["unique.gsm"], ddsoft["soft.unique.gse"],
        ddsoft["soft.unique.gsm"],ddsoft["gsm_soft_and_idat"],
        ddidat[".*idat.gz$"], ddidat[".*idat$"], ddidat[".*hlink.*"],
        ddsoft[ddsoft_lkeys[3]],ddsoft[ddsoft_lkeys[4]],ddsoft[ddsoft_lkeys[5]]
        ]
    }
)

# perc data
val1 = ddidat["eqd.fract.gsm"]
val2 = ddsoft["eqd.fract.gse"]
val3 = ddsoft["eqd.fract.gsm.soft"]
val4 = ddsoft["eqd.fract.gsm.both"]
var_metric = ["GSM IDAT"] * 2
var_metric = var_metric + ["GSE SOFT"] * 2 
var_metric = var_metric + ["GSM SOFT"] * 2 
var_metric = var_metric + ["GSM SOFT & IDAT"] * 2 
var_perc_bp = [val1, 1-val1, val2, 1-val2, val3, 1-val3, val4, 1-val4]
var_perc_df = [val1, val2, val3, val4]
var_type = ["complete", "incomplete"] * 4
dfperc = pd.DataFrame({"type" : [v for v in set(var_metric)],
                        "percent_found" : [100*p for p in var_perc_df]
                        }
                    )
bpdf = pd.DataFrame({"type" : var_metric, 
                    "percent_found" : var_perc_bp,
                    "status" : var_type
                    }
                )


# define browser display
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container(
    html.Div([
        # title
        dbc.Row(dbc.Col(html.Div(
            html.H4('Server dashboard')))),
        # timer
        dbc.Row(dbc.Col(html.Div(id='live-update-text'))),
        # eqperc data
        dbc.Row(dbc.Col(html.Div(dcc.Graph(id='table-eqperc')))),
        dbc.Row(dbc.Col(dcc.Graph(id='table-filecounts'))),
        dbc.Row(dbc.Col(html.Div(dcc.Graph(id='live-barplots-all')))),
        #dbc.Row([
        #        dbc.Col(html.Div(dcc.Graph(id='table-eqperc')), width = 3),
        #        dbc.Col(html.Div(dcc.Graph(id='live-barplots-all')), width = 7)
        #    ]),
        # num files data
        dbc.Row(dbc.Col(dcc.Graph(id='live-scatterplot-all'))),
        dbc.Row(dbc.Col(html.Div(dcc.Interval(id='interval-component',
                    interval=NUM_SEC*1000, # in milliseconds
                    n_intervals=0))))
        ]
    )
)

# TEXT OUTPUTs
# time, text data
@app.callback(Output('live-update-text', 'children'),
              Input('interval-component', 'n_intervals'))
def update_metrics(n):
    time = datetime.datetime.now()
    style = {'padding': '5px', 'fontSize': '16px'}
    mstr = "".join(["Time interval: ", str(NUM_SEC),
        "; elapsed: ", str(n*NUM_SEC), " seconds"])
    return [html.Span(mstr)]

# TABLE OUTPUTS
# percentage found table
@app.callback(Output("table-eqperc", 'figure'),
              Input('interval-component', 'n_intervals'))
def update_table_eqperc(n):
    lval = [c for c in dfperc.columns]
    fig = go.Figure(data=[go.Table(header=dict(values=lval),
                     cells=dict(values=[[d for d in dfperc.type], 
                            [d for d in dfperc.percent_found]]))])
    fig.update_layout(height = 200, width = 300, margin=dict(r=5, l=5, t=20, b=5))
    return fig

# num files table
@app.callback(Output("table-filecounts", 'figure'),
              Input('interval-component', 'n_intervals'))
def update_table_filecounts(n):
    lval = [c for c in dfcount.columns]
    fig = go.Figure(data=[go.Table(header=dict(values=lval),
            cells=dict(values=[
                [d for d in dfcount.file_type], 
                [d for d in dfcount.count_type],
                [d for d in dfcount.counts]
                ]
            )
        )
    ])
    fig.update_layout(height = 300, width=450, 
        margin=dict(r=5, l=5, t=20, b=5))
    return fig

# percent complete, barplots
@app.callback(Output('live-barplots-all', 'figure'),
              Input('interval-component', 'n_intervals'))
def ubarplot_percbp(n):
    """
    """
    figbp = px.bar(bpdf, x="type", y="percent_found", color="status")
    figbp.update_layout(height = 250, width=550, 
        margin=dict(r=5, l=5, t=20, b=5))
    return figbp

# all quantities, overlay
@app.callback(Output('live-scatterplot-all', 'figure'),
              Input('interval-component', 'n_intervals'))
def uscatter_numidat(n):
    """
    """
    # get updated data
    global data_all
    data_all = pd.concat([data_all, 
        pd.DataFrame({"num_idat" : [ddidat['.*idat.gz$']], 
                        "num_gsesoft" : [ddsoft[ddsoft_lkeys[2]]],
                        "num_gsmsoft" : [ddsoft[ddsoft_lkeys[4]]], 
                        "time_seconds" : n*NUM_SEC})])
    fig = go.Figure()
    # numidat
    fig.add_trace(go.Scatter(x=list(data_all.time_seconds), 
        y=list(data_all.num_idat), line_color='#EF553B',
        mode='lines+markers',name='Num. IDATs'))
    # num_gsesoft
    fig.add_trace(go.Scatter(x=list(data_all.time_seconds), 
        y=list(data_all.num_gsesoft), line_color='#19D3F3',
        mode='lines+markers',name='Num. GSE SOFT'))
    # num_gsmsoft
    fig.add_trace(go.Scatter(x=list(data_all.time_seconds), 
        y=list(data_all.num_gsmsoft), line_color='#AB63FA',
        mode='lines+markers',name='Num. GSM SOFT'))
    fig.update_layout(height = 400, width = 650, 
        margin=dict(r=5, l=5, t=20, b=5))
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
