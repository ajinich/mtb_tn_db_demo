import os
from static_vars_functions import *

import dash
import dash_bootstrap_components as dbc
# import dash_core_components as dcc
# import dash_html_components as html
# import dash_table as dt
import dash_daq as daq
from dash import dcc, html
from dash import dash_table as dt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from numpy import inf


# Layout for page analyze datasets
analyze_datasets = html.Div([dbc.Row([html.Label('Pick a dataset')]),
                             dbc.Row([
                                 html.Label("Select by"),
                                 html.Br(),
                                 html.Br(),
                                 dbc.Col([
                                     dcc.Dropdown(id="sel_filter", options=[{"label": x, "value": x} for x in unique_filter_col], value=unique_filter_col[0]),
                                     dcc.Dropdown(id="filter_conditions", clearable=False, value="All"),
                                     dcc.Dropdown(id='sel_dataset',value=unique_expts[0], clearable=False),
                                     dcc.Dropdown(
                                         id='sel_standardized', clearable=False)
                                 ], width=5),
                                 dbc.Col([
                                     daq.Slider(id='log2FC', min=0, max=6, value=1, step=0.5,
                                                size=300,
                                                marks={x: x for x in range(0, 7)}, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "log2FC"})
                                 ], width=3),
                                 dbc.Col([
                                     daq.Slider(id='q-val', min=0, max=1, value=0.05, step=0.05,
                                                marks={
                                                    x / 10: x / 10 for x in range(1, 11)},
                                                size=300, color='#e95420',
                                                handleLabel={"showCurrentValue": True, "label": "q-val"})
                                 ], width=3),
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
                                  for x in unique_Rvs_genes],
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


# Layout for page co-essentiality
co_essentiality = html.Div([
    dbc.Row([html.Label('Pick a gene')]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='sel_gene',
                         options=[{'label': x, 'value': x}
                                  for x in co_genes],
                         placeholder='Select a gene',
                         multi=False,
                         searchable=True, value=co_genes[0]),
            dcc.Dropdown(id='sel_warped_gene',
                         options=[
                             {'label': x, 'value': x} for x in ['Raw Data', 'Warped']],
                         value='Warped', clearable=False,
                         multi=False)
        ]),
        dbc.Col([
            html.Div(id='gene_metadata'),
            html.Br(),
            html.Div([
            html.Div([
                html.Span("Raw Data", style={'color': 'black', 'background-color': '#f0f0f0'}),  # Highlight "Raw Data"
                    ": This term denotes the unprocessed log2 fold-change values derived straight from TnSeq experiments."
                    ]),
            html.Br(),
            html.Div([
            html.Span("Warped Data", style={'color': 'black', 'background-color': '#f0f0f0'}),  # Highlight "Warped Data"
                    ": This refers to the adjusted log2 fold-change values that have been transformed via the Generalized Least Squares method. For a detailed explanation, refer to the study 'Genome-wide co-essentiality analysis in Mycobacterium tuberculosis uncovers an itaconate defense enzyme module' by Jinich et al., published in 2022."
                    ])
])
])
    ], style={'background-color': '#f5f5f5',
              'padding': '30px',
              'border-radius': '25px',
              'border-color': '#dcdcdc',
              'border-width': '2px',
              'border-style': 'solid'}),
    html.Br(),
    html.Br(),
    dbc.Row([
        dbc.Col([
            dcc.Loading(id='loading_correlation', children=dcc.Graph(id='correlation_plot'))
        ])
       ]),
    html.Br(),
    html.Br(),
    dt.DataTable(id='sel_coesen_table',
                 columns=sel_coesen_table_columns,
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
           href='https://github.com/ajinich/mtb_tn_db_demo/blob/master/raw_data_mtb_tn_db.zip?raw=true'),
    # dbc.Button("Download raw data", href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # html.Label('Download raw_data')
])

# Layout for page README
README = html.Div([
    html.P(['MtbTnDB (',
        html.A("www.mtbtndb.app", href="https://www.mtbtndb.app", target="_blank"),
        ") is an online database designed to facilitate the querying and analysis of standardized Mycobacterium tuberculosis (Mtb) Transposon Sequencing (TnSeq) datasets. The platform enables researchers to explore genetic fitness and coessentiality relationships using three primary modes of interaction:"]),
    html.Ul([
        html.Li([
            html.B("Analyze Datasets:", id="adatasets", style={"cursor": "pointer", "color": "blue"}),
            " View and explore TnSeq experimental screens.",
            html.Br(),
            html.Div([
                html.P(["The", html.B("Analyze Datasets"), " tab provides an interactive and information-rich overview of selected TnSeq experimental screens. Each dataset represents a pairwise comparison between an experimental condition and a control."]),
                html.B("Features & Functionality"),
                    html.Ul([
                        html.Li([
                            html.B("About this dataset"),": A summary of the screen including:", 
                            html.Ul([
                                html.Li("Description of the experimental condition."),
                                html.Li("Number of control and experimental replicates."),
                                html.Li("Links to the original publication.")
                                ])]),
                        html.Li([
                            html.B("Volcano Plot & Mutant Interrogation:"),
                            html.Ul([
                                html.Li("Displays differentially represented mutants."),
                                html.Li("Adjustable log2 Fold Change (log2FC) and q-value cutoffs using interactive sliders."),
                                html.Li("Selection of specific genes to highlight in the plot.")
                                ])]),
                        html.Li([
                            html.B("Additional Visualizations"),
                            html.Ul([
                                html.Li([html.B("COG Functional Category Bar Plot"), ": Highlights enriched Clusters of Orthologous Groups (COG) categories."]),
                                html.Li([html.B("Conditional Essentiality vs. Annotation Status Bubble Plot"),": Identifies genes that are conditionally essential but poorly annotated, guiding follow-up studies."])
                                ])]),
                        ]),
                html.B("Example Use Case"),
                html.Br(),
                "A researcher studying Mtb under a stress condition can use the Analyze Datasets tab to:",
                html.Br(),
                html.Ul([
                    html.Li("Select the relevant TnSeq dataset."),
                    html.Li("Examine gene fitness under experimental conditions."),
                    html.Li("Identify genes with significantly altered fitness using volcano and functional category plots.")
                    ]),
                html.Br()
                ], id="adatasets-c", style={"display": "none"})

            ]),
         html.Li([
            html.B("Analyze Genes:", id="header-2", style={"cursor": "pointer", "color": "blue"}),
            " Retrieve gene-specific experimental data.",
            html.Br(),
            html.Div([
                html.P(["The ", html.B("Analyze Genes"), " tab allows users to query a specific gene and view its fitness profile across multiple TnSeq screens."]),
                html.B("Features & Functionality"),
                html.Ul([
                    html.Li([html.B("Gene Information")]),
                    html.Ul([
                        html.Li("Known/probable function."),
                        html.Li("Functional category classification (e.g., Tuberculist annotations)."),
                        html.Li("Essentiality status (Essential/Nonessential, based on mBio criteria).")
                        ]),
                    html.Li([html.B("Experimental Fitness Data")]),
                    html.Ul([
                        html.Li("Log2FC values and q-values for each dataset."),
                        html.Li("Condition-specific fitness profiles."),
                        html.Li([html.B("Track Views"),": Graphical representation of transposon insertion sites within the gene."])
                        ])
                    ]),
                html.B("Example Use Case"),
                html.Br(),
                "A researcher interested in the role of a specific gene under various conditions can:",
                html.Br(),
                html.Ul([
                    html.Li("Search for the gene by name or ID."),
                    html.Li("Review fitness data across multiple screens."),
                    html.Li("Examine essentiality calls and insertion track views to infer potential functional significance.")
                    ]),
                html.Br()
                ], id="content-2", style={"display": "none"})
            ]),
        html.Li([
            html.B("Coessentiality Analysis:", id="header-3", style={"cursor": "pointer", "color": "blue"}),
            " Identify and visualize correlated gene relationships.",
            html.Br(),
            html.Div([
                html.P(["The ", html.B("Coessentiality"), " tab enables users to explore functional relationships between genes by analyzing statistical correlations in TnSeq fitness profiles across multiple experimental conditions."]),
                html.B("Features & Functionality"),
                html.Ul([
                    html.Li([html.B("Gene-Gene Correlation Data"),":"]),
                    html.Ul([
                        html.Li("List of significantly correlated genes with p-values."),
                        html.Li("Known/probable function and functional category of the queried gene.")
                        ]),
                    html.Li(html.B("Interactive Visualization:")),
                    html.Ul([
                        html.Li("Graphical representation of coessential gene networks."),
                        html.Li("Displays both first-degree and second-degree relationships, allowing exploration of broader genetic interactions.")
                        ]),
                    ]),
                html.B("Example Use Case"),
                html.Br(),
                "A researcher hypothesizing that a gene interacts with another can:",
                html.Br(),
                html.Ul([
                    html.Li("Input the gene name or ID."),
                    html.Li("Review correlated genes and their statistical significance."),
                    html.Li("Visualize the gene network to infer functional dependencies.")]),
                    html.Br()                
                ], id="content-3", style={"display": "none"})

            ]),
        html.Li([
            html.B("Data Accessibility & Download Options", id="header-4", style={"cursor": "pointer", "color": "blue"}),
            html.Div(["Users can download individual TnSeq screens from the ", html.B("Analyze Datasets"), " tab. The ", html.B("About"), " tab provides details on data sources and updates."], id="content-4", style={"display": "none"})
            ])
        ]),
    html.P('The modular design of MtbTnDB supports the addition of future datasets, including conditional essentiality screens using CRISPRi.'),
    html.B("Conclusion"),
    html.P("MtbTnDB provides a powerful and intuitive interface for exploring Mtb TnSeq data. By integrating interactive visualizations, standardized fitness data, and coessentiality analysis, the platform enables researchers to uncover gene functions, essentiality patterns, and functional relationships in Mtb. Future updates will expand dataset coverage, incorporating new methodologies such as CRISPRi-based screens.")
])
