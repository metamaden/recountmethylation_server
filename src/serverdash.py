#!/usr/bin/env python3

# Author: Sean Maden
# 
# Dashboard script for recountmethylation server.
#
#

import datetime, dash, plotly, os, sys, re
import pandas as pd
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
from report import idats_report, soft_report

# GLOBALS
# define globals
NUM_SEC = 5
# new query
didat = idats_report()
dsoft = soft_report(didat = didat)
# main data for line plots by time
data_all = pd.DataFrame({"num_idat" : [0], "num_gsesoft" : [0],
    "num_gsmsoft" : [0], "time" : [0]})
# file counts data
dfcount = pd.DataFrame({"file_type" : ["IDATs (all)", "IDATs (only expanded)",
                                    "GSE SOFT (all)", "GSE SOFT (only expanded)",
                                    "GSM SOFT (all)", "GSM SOFT (with paired IDAT)",
                                    "GSM SOFT (with JSON)", "GSM JSON (all)"],
                        "num_files" : [didat["total_files"], 
                                        didat["total_expanded_files"],
                                        dsoft["soft_gse_total_files"],
                                        dsoft["soft_gse_total_expanded_files"],
                                        dsoft["soft_gsm_total_files"],
                                        dsoft["soft_num_gsm_ipaircomp"],
                                        dsoft["gsm_soft_with_json"],
                                        dsoft["json_gsm_total_files"]
                                    ]
                        }
                    )

# perc data
val1 = round(didat["expanded_fract_of_eqidats"], 3)
val2 = round(dsoft["soft_gse_all_fract_target"], 3)
val3 = round(dsoft["soft_gsm_all_fract_target"], 3)
val4 = round(dsoft["json_gsm_fract_target"], 3)
var_metric = ["IDAT"] * 2
var_metric = var_metric + ["GSE SOFT"] * 2 
var_metric = var_metric + ["GSM SOFT"] * 2 
var_metric = var_metric + ["GSM JSON"] * 2 
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
            html.H4('Recountmethylation server dashboard')))),
        # timer
        dbc.Row(dbc.Col(html.Div(id='live-update-text'))),
        # num files data
        dbc.Row([
                dbc.Col(dcc.Graph(id='table-filecounts'), width = 5),
                dbc.Col(dcc.Graph(id='live-scatterplot-all'), width = 6)
            ]),
        # eqperc data
        dbc.Row([
                dbc.Col(html.Div(dcc.Graph(id='table-eqperc')), width = 5),
                dbc.Col(html.Div(dcc.Graph(id='live-barplots-all')), width = 6)
            ]),
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
# num files table
@app.callback(Output("table-filecounts", 'figure'),
              Input('interval-component', 'n_intervals'))
def update_table_filecounts(n):
    lval = [c for c in dfcount.columns]
    fig = go.Figure(data=[go.Table(header=dict(values=lval),
            cells=dict(values=[[d for d in dfcount.file_type], 
                                [d for d in dfcount.num_files]]
                                )
            )
    ])
    fig.update_layout(height = 250, width=450, 
        margin=dict(r=5, l=5, t=20, b=5))
    return fig

# percentage found table
@app.callback(Output("table-eqperc", 'figure'),
              Input('interval-component', 'n_intervals'))
def update_table_eqperc(n):
    lval = [c for c in dfperc.columns]
    fig = go.Figure(data=[go.Table(header=dict(values=lval),
                     cells=dict(values=[[d for d in dfperc.type], 
                            [d for d in dfperc.percent_found]]))])
    fig.update_layout(height = 300, margin=dict(r=5, l=5, t=20, b=5))
    return fig

# PLOT OUTPUTS
# percent complete, barplots
@app.callback(Output('live-barplots-all', 'figure'),
              Input('interval-component', 'n_intervals'))
def ubarplot_percbp(n):
    """
    """
    figbp = px.bar(bpdf, x="type", y="percent_found", color="status")
    figbp.update_layout(height = 250, width=450, margin=dict(r=5, l=5, t=20, b=5))
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
        pd.DataFrame({"num_idat" : [didat["total_files"]], 
                        "num_gsesoft" : [dsoft["soft_gse_total_files"]],
                        "num_gsmsoft" : [dsoft["soft_gsm_total_files"]], 
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
    fig.update_layout(height = 200, margin=dict(r=5, l=5, t=20, b=5))
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
