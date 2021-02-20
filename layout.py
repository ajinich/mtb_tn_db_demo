import os
from static_vars_functions import *

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_daq as daq
import dash_html_components as html
import dash_table as dt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from numpy import inf

# Layout for page analyze datasets
analyze_datasets = html.Div([dbc.Row([html.Label('Pick a dataset')]),
                             dbc.Row([
                                 html.Br(),
                                 html.Br(),
                                 dbc.Col([
                                     dcc.Dropdown(id='sel_dataset',
                                                  options=[{'label': x, 'value': x}
                                                           for x in unique_expts],
                                                  value=unique_expts[0], clearable=False),
                                     dcc.Dropdown(
                                         id='sel_standardized', clearable=False)
                                 ], width=4),
                                 dbc.Col([
                                     daq.Slider(id='log2FC', min=0, max=6, value=1, step=0.5,
                                                size=300,
                                                marks={x: x for x in range(0, 7)}, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "log2FC"})
                                 ], width=4),
                                 dbc.Col([
                                     daq.Slider(id='q-val', min=0, max=1, value=0.05, step=0.05,
                                                marks={
                                                    x / 10: x / 10 for x in range(1, 11)},
                                                size=300, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "q-val"})
                                 ], width=4),
                                 html.Br(),
                                 html.A('Download this dataset', id='download_dataset', download="", href="",
                                        target="_blank"),

                             ], align='center',
                                 style={'background-color': '#f5f5f5',
                                        'padding': '30px',
                                        'border-radius': '25px',
                                        'border-color': '#dcdcdc',
                                        'border-width': '2px',
                                        'border-style': 'solid'}),
                             html.Br(),
                             html.Br(),
                             dbc.Row([
                                 dbc.Col([
                                     html.Div([
                                         html.Label('About this dataset')
                                     ], style={'textAlign': 'center', 'display': 'block'}),
                                 ], align='center', width=3),
                                 dbc.Col([
                                     html.Div([
                                         html.Label('Volcano plot')
                                     ], style={'textAlign': 'center', 'display': 'block'}),
                                     html.Div(id='num_significant', style={
                                         'textAlign': 'center', 'display': 'block'}),
                                 ], align='center', width=5),
                                 dbc.Col([
                                     html.Div([
                                         html.Label('Gene List')
                                     ], style={'textAlign': 'center', 'display': 'block'}),
                                 ], align='center', width=3),
                             ]),
                             dbc.Row([
                                 dbc.Col([
                                     html.Div(id='dataset_metadata')], width=3,
                                     style={'background-color': '#f5f5f5',
                                            'padding': '30px',
                                            'border-radius': '25px',
                                            'border-color': '#dcdcdc',
                                            'border-width': '2px',
                                            'border-style': 'solid'}),
                                 dbc.Col([dcc.Loading(id='loading_volcano', children=dcc.Graph(id='volcano')),
                                          ],
                                         width=5, align='center'),
                                 dbc.Col([dcc.Loading(id='loading_dataset_table',
                                                      children=dt.DataTable(id='sel_dataset_table',
                                                                            columns=[{"name": i, "id": i} for i in [
                                                                                'Rv_ID', 'gene_name', 'log2FC',
                                                                                'q-val']],
                                                                            sort_action='native',
                                                                            row_selectable='multi',
                                                                            selected_rows=[],
                                                                            page_action='native',
                                                                            page_size=15,
                                                                            page_current=0,
                                                                            style_header={'color': '#e95420',
                                                                                          'font-weight': 'bold',
                                                                                          'text-align': 'center'},
                                                                            style_cell_conditional=[
                                                                                {'if': {'column_id': 'q-val'},
                                                                                 'width': '30%'}
                                                                            ],
                                                                            style_data_conditional=[
                                                                                {'if': {'row_index': 'odd'},
                                                                                 'backgroundColor': 'rgb(248,248,248)'}
                                                                            ],
                                                                            style_cell={
                                                                                'font-family': 'ubuntu',
                                                                                'font-size': 14,
                                                                                'height': '10px',
                                                                                'textOverFlow': 'ellipsis',
                                                                                'text-align': 'center',
                                                                                'overflow': 'hidden'
                                                                            },
                                                                            #    style_as_list_view=True,
                                                                            ))
                                          ], width=4, align='center')
                             ]),
                             html.Br(),
                             html.Br(),
                             html.Br(),
                             dbc.Row([
                                 dbc.Col([
                                     html.Div([
                                         html.Label('COG Categories')
                                     ], style={'textAlign': 'center', 'display': 'block'}),
                                 ], align='center', width=6),
                                 dbc.Col([
                                     html.Div([
                                         html.Label('Essentiality plot')
                                     ], style={'textAlign': 'center', 'display': 'block'}),
                                 ], align='center', width=6),
                             ]),
                             dbc.Row([
                                 dbc.Col([
                                     dcc.Dropdown(id='sel_cog',
                                                  options=[{'label': x, 'value': x} for x in
                                                           ['Under-represented', 'Over-represented']],
                                                  value='Under-represented'),
                                     dcc.Graph(id='cog')

                                 ], width=6, align='center'),
                                 dbc.Col([dcc.Loading(id='loading_bubble', children=dcc.Graph(id='bubble_plot'))
                                          ], width=6, align='center')
                             ], justify='center')
                             ])

# Layout for page analyze genes
analyze_genes = html.Div([
    dbc.Row([html.Label('Pick a gene')]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='sel_gene',
                         options=[{'label': x, 'value': x}
                                  for x in unique_genes + unique_Rvs],
                         placeholder='Select a gene',
                         multi=False,
                         searchable=True),
            dcc.Dropdown(id='sel_standardized_gene_table',
                         options=[
                             {'label': x, 'value': x} for x in ['Standardized', 'Original']],
                         value='Standardized', clearable=False,
                         multi=False)
        ]),
        dbc.Col([
            html.Div(id='gene_metadata')])
    ], style={'background-color': '#f5f5f5',
              'padding': '30px',
              'border-radius': '25px',
              'border-color': '#dcdcdc',
              'border-width': '2px',
              'border-style': 'solid'}),
    html.Br(),
    html.Br(),
    dt.DataTable(id='sel_gene_table',
                 columns=sel_gene_table_columns,
                 sort_action='native',
                 page_action='native',
                 page_size=15,
                 page_current=0,
                 style_cell={
                     'font-family': 'ubuntu',
                     'font-size': 14,
                     'whiteSpace': 'normal',
                     'height': 'auto',
                     'maxWidth': '120px',
                     # 'textOverFlow': 'ellipsis',
                     'textAlign': 'center',
                     #  'overflow': 'hidden'
                 },
                 style_header={'color': '#e95420', 'font-weight': 'bold',
                               'textAlign': 'center'},
                 style_data_conditional=[
                     {'if': {'row_index': 'odd'},
                      'backgroundColor': 'rgb(248,248,248)'},
                     {'if': {'column_id': 'Expt'}, 'maxWidth': '300px'},
                     {'if': {'column_id': 'meaning'}, 'width': '250px'},
                     {'if': {'column_id': 'log2FC'}, 'width': '60px'},
                     {'if': {'column_id': 'q-val'}, 'width': '60px'},
                     {'if': {'column_id': 'gene_name'}, 'width': '100px'},
                     {'if': {'column_id': 'num replicates control'}, 'width': '50px'},
                     {'if': {'column_id': 'num replicates experimental'},
                         'width': '50px'},
                     {'if': {'column_id': 'Expt'}, 'whiteSpace': 'pre-line'},
                     {'if': {'column_id': 'paper'}, 'width': '300px'},
                 ],
                 )

])

# Layout for page About
about = html.Div([
    html.Span('TnSeq has been used extensively in '),
    html.Span('M. tuberculosis', style={'font-style': 'italic'}),
    html.Span(' genetic research and identification of gene essentiality (TnSeq) profiles is important for predicting gene function. However, these profiles are buried across dozens of research papers within supplementary materials which makes querying them cumbersome. The MtbTnDB solves this problem by building a central repository of TnSeq screens performed in '),
    html.Span('M. tuberculosis', style={'font-style': 'italic'}),
    html.Span(
        ', and allows users easy access to data through an interactive web-app.'),
    html.Br(),
    html.Br(),
    html.H5('Contact'),
    html.Span("For bug reports and data submissions, contact "),
    html.A('Adrian Jinich',
           href="mailto:adj2010@med.cornell.edu", target='_blank'),
    html.Br(),
    html.Br(),
    html.H5('Raw data'),
    html.Span('Raw data is available '),
    html.A('here',
           href='https://www.google.com'),
    # dbc.Button("Download raw data", href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # html.Label('Download raw_data')
])
