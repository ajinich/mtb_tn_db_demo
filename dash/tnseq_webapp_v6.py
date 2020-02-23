import os
import urllib.parse
from io import StringIO
from random import random

import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_daq as daq
import dash_html_components as html
import dash_table as dt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import seaborn as sns
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from numpy import inf

external_stylesheets = [dbc.themes.UNITED]
path_data = '../data/'
main_data = pd.read_csv(os.path.join(path_data, 'data_dash_meta.csv'))
path_annotation = '../data/annotations/'
cogs_df = pd.read_csv(os.path.join(path_annotation, 'all_cogs.csv'))
cogs_desc = pd.read_csv(os.path.join(
    path_annotation, 'cog_names.csv'), header=0, index_col=0, squeeze=True)

unique_expts = list(main_data['Expt'].unique())
unique_Rvs = sorted(list(main_data['Rv_ID'].unique()))
unique_genes = sorted(list(main_data['gene_name'].unique()))
main_data['id'] = main_data['Rv_ID']
main_data.set_index('id', inplace=True, drop=False)
# print("Main", main_data.head())
df_uk = pd.read_csv(os.path.join(
    path_annotation, 'unknown_essentials/unknown_ALL_levels_essential_scores.csv'))
df_uk = df_uk[['Rv_ID', 'gene_name', 'UK_score_4']]


def discretize_q_values(row):
    q_val = row['q-val']
    if q_val < 0.01:
        q_val_d = 3
    elif q_val < 0.05:
        q_val_d = 2
    else:
        q_val_d = 1
    return q_val_d


def unknown_essential_xy(TnSeq_screen, rand_param=0.6):
    # Grab data for a single TnSeq screen
    selected_data = main_data[main_data['Expt'] == TnSeq_screen]
    selected_data['q_val_D'] = selected_data.apply(discretize_q_values, 1)

    # Merge with unknowns:
    df_vis = selected_data.merge(df_uk, on=['Rv_ID'], how='inner')

    # Get x-y datasets:
    rv_ids = df_vis.Rv_ID.values
    uk_list = np.array(df_vis.UK_score_4)
    q_list = np.array(df_vis.q_val_D)

    # randomize:
    uk_rd = np.array([uk + rand_param * random() -
                      rand_param / 2 for uk in uk_list])
    q_rd = np.array([q + rand_param * random() -
                     rand_param / 2 for q in q_list])

    # color the unknown-essentials differently:
    color_list = np.array(['#585858']*df_vis.shape[0])
    color_list[df_vis[(df_vis.q_val_D == 3) &
                      (df_vis.UK_score_4 == 4)].index] = '#2b7bba'
    return uk_rd, q_rd, color_list, rv_ids


app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

navbar = dbc.NavbarSimple([
    dbc.NavItem(dbc.NavLink('Analyze datasets',
                            active=True, href=app.get_relative_path('/analyze_datasets'))),
    dbc.NavItem(dbc.NavLink('Analyze genes',
                            href=app.get_relative_path('/analyze_genes'))),
    dbc.NavItem(dbc.NavLink('About', active=True,
                            href=app.get_relative_path('/about')))
], brand="Mtb Tn-seq database", color='primary', light=True)

analyze_datasets = html.Div([dbc.Row([html.Label('Pick a dataset')]),
                             dbc.Row([
                                 html.Br(),
                                 html.Br(),
                                 dbc.Col([
                                     dcc.Dropdown(id='sel_dataset',
                                                  options=[{'label': x, 'value': x}
                                                           for x in unique_expts],
                                                  value=unique_expts[0])
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
        dbc.Col([dcc.Graph(id='volcano'),
                 ],
                width=5, align='center'),
        dbc.Col([dt.DataTable(id='sel_dataset_table',
                              columns=[{"name": i, "id": i} for i in [
                                  'Rv_ID', 'gene_name', 'log2FC', 'q-val']],
                              sort_action='native',
                              row_selectable='multi',
                              selected_rows=[],
                              page_action='native',
                              page_size=15,
                              page_current=0,
                              style_header={'color': '#e95420', 'font-weight': 'bold',
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
                              style_as_list_view=True,
                              )
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
            dcc.Dropdown(id='Sel_cog',
                         options=[{'label': x, 'value': x} for x in
                                  ['Under-represented', 'Over-represented']],
                         value='Under-represented'),
            dcc.Graph(id='cog')

        ], width=6, align='center'),
        dbc.Col([
            dcc.Graph(id='bubble_plot')
        ], width=6, align='center')
    ], justify='center')
])

analyze_genes = html.Div([
    dbc.Row([html.Label('Pick a gene')]),
    dbc.Row([
        dbc.Col([
            dcc.Dropdown(id='sel_gene', options=[{'label': x, 'value': x} for x in unique_genes+unique_Rvs],
                         placeholder='Select a gene', multi=False, searchable=True)]),
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
                 columns=[{"name": i,
                           "id": i,
                           "presentation": 'markdown'} for i in [
                     'Rv_ID', 'gene_name', 'Expt', 'log2FC', 'q-val', '# control reps', '# experimental reps']],
                 sort_action='native',
                 page_action='native',
                 #  filter_action='native',
                 page_size=15,
                 page_current=0,
                 style_header={'color': '#e95420', 'font-weight': 'bold',
                               'text-align': 'center'},
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
                 style_as_list_view=True)
])

about = html.Div([
    html.Label('Desc and credits'),
    html.Br(),
    html.A('Download raw data',
           href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # dbc.Button("Download raw data", href='https://www.dropbox.com/s/ktx859tq73i8y9m/ORF_details_final.csv?dl=1'),
    # html.Label('Download raw_data')

])

app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        navbar,
        dbc.Container(id="content", style={"padding": "20px"}),
    ])

app.config.suppress_callback_exceptions = True
app.scripts.config.serve_locally = True


@app.callback(Output('loading_dataset_table', "children"))
@app.callback(Output("loading_volcano", "children"))
@app.callback(Output("content", "children"),
              [Input("url", "pathname")])
def display_content(path):
    page_name = app.strip_relative_path(path)
    if page_name == 'analyze_datasets':
        return analyze_datasets
    if page_name == "analyze_genes":
        return analyze_genes
    if page_name == 'about':
        return about


@app.callback([
    Output('download_dataset', 'href'),
    Output('download_dataset', 'download')
],
    [Input('sel_dataset', 'value')])
def update_download_dataset(sel_dataset):
    dff = main_data[main_data['Expt'] == sel_dataset]
    csv_string = dff.to_csv(encoding='utf-8')
    csv_string = "data:text/csv;charset=utf-8," + \
        urllib.parse.quote(csv_string)
    download_string = sel_dataset + '.csv'
    return csv_string, download_string


@app.callback(
    Output('dataset_metadata', 'children'),
    [Input('sel_dataset', 'value')])
def print_dataset_metadata(sel_dataset):
    dff = main_data[main_data['Expt'] == sel_dataset]
    text = [
        html.Strong('Summary'),
        html.Span(': ' + dff['meaning'][0]),
        html.Br(),
        html.Br(),
        html.Strong('Original Publication:'),
        html.Br(),
        html.A(dff['paper_title'][0],
               href=dff['paper_URL'][0]),
        html.Br(),
        html.Br(),
        html.Strong('No of control replicates'),
        html.Span(': ' + str(dff['# control reps'][0])),
        html.Br(),
        html.Br(),
        html.Strong('No of experimental replicates'),
        html.Span(': ' + str(dff['# experimental reps'][0]))
    ]
    return text


@app.callback(
    Output('volcano', 'figure'),
    [Input('sel_dataset', 'value'),
     Input('log2FC', 'value'),
     Input('q-val', 'value'),
     Input('sel_dataset_table', "derived_virtual_selected_row_ids")])
def update_volcano(sel_dataset, log2FC, qval, selected_row_ids):
    dff = main_data[main_data['Expt'] == sel_dataset]
    if selected_row_ids is None:
        selected_row_ids = []

    # make qval ticks, replacing the np.nans with inf
    max_log_qval = np.unique(-np.log10(dff['q-val']))[-2]
    inf_repl = np.ceil(max_log_qval) + 1
    dff['qval_plotting'] = -np.log10(dff['q-val'])
    dff['qval_plotting'].replace(np.inf, inf_repl, inplace=True)
    tickvals = list(np.arange(0, inf_repl + 0.5, 0.5))
    ticklab = tickvals.copy()
    ticklab[-1] = 'Inf'

    for_x_ticks = dff['log2FC']
    # TODO: WHAT is this?
    for_x_ticks = for_x_ticks.replace([np.inf, -np.inf], np.nan)
    for_x_ticks = for_x_ticks.dropna()
    tickvals_x = list(np.arange(int(for_x_ticks.min() - 1),
                                int(for_x_ticks.max() + 1), 1))
    ticklab_x = tickvals_x.copy()

    # split data into tickeed,unticked_sig, unticked_non_sig
    ticked = dff['id'].isin(selected_row_ids)
    ticked_data = dff[ticked]
    unticked_data = dff[~ticked]
    generated_filter = (unticked_data['q-val'] <= qval) & (
        (unticked_data['log2FC'] <= (-log2FC)) | (unticked_data['log2FC'] >= log2FC))
    sig_data = unticked_data[generated_filter]
    non_sig_data = unticked_data[~generated_filter]

    traces = []
    traces.append(go.Scatter(x=sig_data['log2FC'],
                             y=sig_data['qval_plotting'],
                             text=sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Outside cutoff',
                             marker={'opacity': 0.6, 'size': 10,
                                       'color': 'orangered'},
                             showlegend=False,
                             ))
    traces.append(go.Scatter(x=non_sig_data['log2FC'],
                             y=non_sig_data['qval_plotting'],
                             text=non_sig_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers',
                             name='Pass cutoff',
                             marker={'opacity': 0.6,
                                       'size': 10,
                                       'color': 'grey'},
                             showlegend=False
                             ))
    traces.append(go.Scatter(x=ticked_data['log2FC'],
                             y=ticked_data['qval_plotting'],
                             text=ticked_data['Rv_ID'],
                             hoverinfo='text',
                             mode='markers+text',
                             textposition='bottom center',
                             name='T',
                             marker={'opacity': 0.6,
                                       'size': 10,
                                       'color': 'green'},
                             showlegend=False
                             ))
    return {'data': traces,
            'layout': go.Layout(
                autosize=False,
                margin={
                    'l': 45,
                    'r': 15,
                    'pad': 0,
                    't': 30,
                    'b': 90, },
                xaxis={'title': 'log2FC', 'ticktext': ticklab_x,
                       'tickvals': tickvals_x},
                yaxis={
                    'title': '-log10(q-val)', 'ticktext': ticklab, 'tickvals': tickvals},
                hovermode='closest'
            )}


@app.callback(
    Output('sel_dataset_table', 'data'),
    [Input('sel_dataset', 'value')])
def update_dataset_table(sel_dataset):
    dff = main_data[main_data['Expt'] == sel_dataset]
    dff = dff[[
        'Rv_ID', 'gene_name', 'log2FC', 'q-val', 'id']]
    dff['q-val'] = np.round(dff['q-val'], 2)
    dff['log2FC'] = np.round(dff['log2FC'], 2)
    dff = dff.sort_values(by='log2FC')
    return dff.to_dict('records')


@app.callback(Output('cog', 'figure'),
              [Input('sel_dataset', 'value'),
               Input('Sel_cog', 'value'),
               Input('log2FC', 'value'),
               Input('q-val', 'value')])
def update_cog(sel_dataset, sel_cog, log2FC, qval):
    # filter data based on user inputs
    selected_data = main_data[main_data['Expt'] == sel_dataset]
    if sel_cog == 'Under-represented':
        sel_subset_filter = (
            selected_data['q-val'] <= qval) & (selected_data['log2FC'] <= -log2FC)
        colorscale = 'Cividis'
    else:
        sel_subset_filter = (
            selected_data['q-val'] <= qval) & (selected_data['log2FC'] >= -log2FC)
        colorscale = 'Viridis'
    sel_subset = selected_data[sel_subset_filter]
    # calculate genome freq
    cog_total_freq = cogs_df['COG'].value_counts(normalize=True)
    # calculate subset freq
    sel_cogs = cogs_df[cogs_df['X.Orf'].isin(sel_subset['Rv_ID'])]
    sel_cogs_freq = sel_cogs['COG'].value_counts(normalize=True)
    # calculate enrichment
    normalized_cogs = sel_cogs_freq / cog_total_freq
    # format
    normalized_cogs = normalized_cogs[~normalized_cogs.isnull()]
    normalized_cogs = normalized_cogs.sort_values()
    cog_names = cogs_desc.loc[list(normalized_cogs.index)]
    cog_names = list(cog_names.values)

    bar_data = [go.Bar(y=list(normalized_cogs.index), x=list(normalized_cogs.values),
                       orientation='h',
                       text=cog_names,
                       hoverinfo='text',
                       marker={'color': list(normalized_cogs.values), 'colorscale': colorscale})]
    return {'data': bar_data,
            'layout': go.Layout(
                margin={
                    'l': 50,
                    'r': 10,
                    'pad': 3,
                    't': 30,
                    'b': 90, },
                # paper_bgcolor='rgba(0,0,0,0)',
                # plot_bgcolor = 'rgba(0,0,0,0)',
                xaxis={'title': 'Normalized to genomic frequency'},
                hovermode='closest',
                shapes=[{'type': 'line', 'x0': 1, 'y0': 0, 'x1': 1, 'y1': len(normalized_cogs),
                         'line': {'color': 'grey', 'width': 1, 'dash': 'dot'}}])
            }


@app.callback(
    Output('bubble_plot', 'figure'),
    [Input('sel_dataset', 'value')])
def update_bubble(sel_dataset):
    uk_rd, q_rd, color_list, rv_ids = unknown_essential_xy(
        sel_dataset)
    print(color_list)
    return {
        'data': [
            go.Scatter(
                x=uk_rd,
                y=q_rd,
                text=rv_ids,
                mode='markers',
                opacity=1,
                hoverinfo='text',
                marker={
                    'size': 15,
                    'line': {'width': 1.5, 'color': 'black'},
                    'color': color_list
                }
            )
        ],
        'layout': go.Layout(
            autosize=True,
            # width=800,
            # height=500,
            xaxis=go.layout.XAxis(
                tickmode='array',
                tickvals=[0, 1, 2, 3, 4],
                ticktext=['most well<br>characterized', '', '', '',
                          'least well<br>characterized'],
                tickfont=dict(size=14),
                title='Annotation'

            ),
            yaxis=go.layout.YAxis(
                tickmode='array',
                tickvals=[1, 2, 3],
                ticktext=['non-essential', 'q-val < 0.05', 'q-val < 0.01'],
                tickangle=270,
                tickfont=dict(size=14),
                title='Essentiality'
            ),
            margin={'l': 30, 'b': 100, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest'
        )
    }


@app.callback(
    Output('sel_gene_table', 'data'),
    [Input('sel_gene', 'value')])
def update_genes_table(selected_gene):
    if selected_gene in unique_Rvs:
        dff = main_data[main_data['Rv_ID'] == selected_gene]
    elif selected_gene in unique_genes:
        dff = main_data[main_data['gene_name'] == selected_gene]
    else:
        raise PreventUpdate
    dff = dff[['Rv_ID', 'gene_name', 'Expt', 'log2FC', 'q-val',
               '# control reps', '# experimental reps']]
    dff['q-val'] = np.round(dff['q-val'], 2)
    dff['log2FC'] = np.round(dff['log2FC'], 2)
    dff = dff.sort_values(by='log2FC')
    return dff.to_dict('records')


@app.callback(
    Output('gene_metadata', 'children'),
    [Input('sel_gene', 'value')])
def print_gene_metadata(sel_gene):
    if sel_gene in unique_Rvs:
        sel_details = main_data[main_data['Rv_ID'] == sel_gene]
    elif sel_gene in unique_genes:
        sel_details = main_data[main_data['gene_name'] == sel_gene]
    else:
        return "gene not found"
    return list(sel_details['Description'])[0]


if __name__ == '__main__':
    app.run_server(debug=True)
