import os
import urllib.parse
from io import StringIO

import dash
import dash_bootstrap_components as dbc
# import dash_core_components as dcc
# import dash_html_components as html
# import dash_table as dt
from dash import dcc, html
from dash import dash_table as dt
import dash_daq as daq
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from numpy import inf
import itertools
from plotly.subplots import make_subplots
from random import sample

# TODO: reformat analyze genes to remove underscores in header, add ellipsis to overflow, esp for paper
#####
# SECTION 1: Read in data, create static variables
#####
# from static_vars_functions import *

#####
# SECTION 2: LAYOUT
from layout import *

#####
# initialize app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
server = app.server

# create navbar
navbar = dbc.NavbarSimple([
    dbc.NavItem(dbc.NavLink('Analyze datasets',
                            active=True, href=app.get_relative_path('/analyze_datasets'))),
    dbc.NavItem(dbc.NavLink('Analyze genes',
                            href=app.get_relative_path('/analyze_genes'))),
    dbc.NavItem(dbc.NavLink('Co-essentiality',
                            href=app.get_relative_path('/co-essentiality'))),
    dbc.NavItem(dbc.NavLink('About', active=True,
                            href=app.get_relative_path('/about')))
], brand="MtbTnDB", color='primary', light=True)

# app.layout dynamically takes in different content based on the path. See next callback
app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        navbar,
        dbc.Container(id="content", style={"padding": "20px"}),
    ])

app.config.suppress_callback_exceptions = True
app.scripts.config.serve_locally = True

#####
# SECTION 3: CALLBACKS
#####


@ app.callback(Output("content", "children"),
               [Input("url", "pathname")])
def display_content(path):
    """
    Takes in path from the URL and returns layout for one of four pages
    """
    page_name = app.strip_relative_path(path)
    if page_name == 'analyze_datasets':
        return analyze_datasets
    if page_name == "analyze_genes":
        return analyze_genes
    if page_name == "co-essentiality":
        return co_essentiality
    if page_name == 'about':
        return about


@ app.callback(
    [Output('sel_standardized', 'options'),
     Output('sel_standardized', 'value')],
    [Input('sel_dataset', 'value')])
def update_standardized_dropdown(sel_dataset):
    """Take in the dataset selected and update the sel_standardized dropdown

    Args:
        sel_dataset (str): Selected dataset eg: griffin_glycerol_vs_mbio_H37Rv

    Returns:
        list: List of options for sel_standardized
        str: Default value to display from the options
    """
    # is it a std_dataset?
    if sel_dataset in dict_std_to_si:
        # does a corresponding original (aka si) dataset exist?
        if pd.isna(dict_std_to_si[sel_dataset]):
            return [{'label': x, 'value': x} for x in ['Standardized']], 'Standardized'
        else:
            return [{'label': x, 'value': x} for x in ['Standardized', 'Original']], 'Standardized'
    # is it an original dataset?
    else:
        # does a corresponding std_dataset exist?
        if pd.isna(dict_si_to_std[sel_dataset]):
            return [{'label': x, 'value': x} for x in ['Original']], 'Original'
        else:
            return [{'label': x, 'value': x} for x in ['Standardized', 'Original']], 'Original'


@ app.callback([
    Output('download_dataset', 'href'),
    Output('download_dataset', 'download'),
    Output('num_significant', 'children'),
],
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value'),
     Input('log2FC', 'value'),
     Input('q-val', 'value'),
     ])
def update_multiple_outputs_analyze_datasets(sel_dataset, sel_standardized, log2FC, qval):
    """Using user inputs, update download dataset and num significant

    Args:
        sel_dataset (str): User selected dataset 
        sel_standardized (str): User selected standardized/original
        log2FC (float): User selected log2FC cutoff
        qval (flaot): User selected qval cutoff

    Returns:
        str: href for download
        str: file name for download
        list: list of text for number of significant. 
              ' ' is returned if not enough genes in experiment for this to be meaningful
    """
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    num_significant_text = update_num_significant(dff, log2FC, qval)
    if sel_standardized == 'Original':
        if dict_plot_si[dataset_name] == 'No':
            num_significant_text = ' '
    csv_string, download_string = update_download_dataset(dff, dataset_name)

    return csv_string, download_string, num_significant_text



@app.callback([Output('volcano', 'figure'),
               Output('volcano', 'config')],
              [Input('sel_dataset', 'value'),
               Input('sel_standardized', 'value'),
               Input('log2FC', 'value'),
               Input('q-val', 'value'),
               Input('sel_dataset_table', "derived_virtual_selected_row_ids")
               ])
def update_volcano(sel_dataset, sel_standardized, log2FC, qval, selected_row_ids):
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    config = {
        'modeBarButtonsToRemove': plotly_buttons_remove,
        'toImageButtonOptions': {
            'height': 700,
            'width': 700,
            'scale': 5,
            'filename': f'{dataset_name}_volcano.png'
        }
    }
    if sel_standardized == 'Original':
        # Is there enough data for a meaningful plot?
        if dict_plot_si[dataset_name] == 'No':
            return (empty_plot('Not enough datapoints' + '\n' + 'for a meaningful plot'), config)
    # weird plotly requirement
    if selected_row_ids is None:
        selected_row_ids = []
    # TODO: Figure out NAs
    dff = dff.dropna(axis=0, subset=['log2FC', 'q-val'])
    # make qval ticks, replacing the np.nans with inf
    # what is current second highest max log10 transformed qval?
    # note that max will be inf
    max_log_qval = np.unique(-np.log10(dff['q-val']))[-2]
    # create new column in dff with qval for plotting, replace inf values
    inf_repl = np.ceil(max_log_qval) + 1
    dff['qval_plotting'] = -np.log10(dff['q-val'])
    dff['qval_plotting'].replace(np.inf, inf_repl, inplace=True)

    # create x and y tick vals and labels
    tickvals = list(np.arange(0, inf_repl + 0.5, 0.5))
    ticklab = tickvals.copy()
    ticklab[-1] = 'Inf'
    for_x_ticks = dff['log2FC']
    # TODO: WHAT is this? - commented for now
    # for_x_ticks = for_x_ticks.replace([np.inf, -np.inf], np.nan)
    # for_x_ticks = for_x_ticks.dropna()
    tickvals_x = list(np.arange(int(for_x_ticks.min() - 1),
                                int(for_x_ticks.max() + 1), 1))
    ticklab_x = tickvals_x.copy()

    # split data into selected (ie ticked), unselected_sig, unselected_non_sig
    ticked = dff['id'].isin(selected_row_ids)
    ticked_data = dff[ticked]
    unticked_data = dff[~ticked]
    generated_filter = (unticked_data['q-val'] <= qval) & (
        (unticked_data['log2FC'] <= (-log2FC)) | (unticked_data['log2FC'] >= log2FC))
    sig_data = unticked_data[generated_filter]
    non_sig_data = unticked_data[~generated_filter]

    # make traces for each kind of data
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
    # return dict of plotting and config for plotly
    return ({'data': traces,
             'layout': go.Layout(
                 autosize=False,
                 margin={'l': 45, 'r': 15, 'pad': 0, 't': 30, 'b': 90},
                 xaxis={'title': 'log2FC', 'ticktext': ticklab_x,
                        'tickvals': tickvals_x, 'fixedrange': True},
                 yaxis={'title': '-log10(q-val)', 'ticktext': ticklab,
                        'tickvals': tickvals, 'fixedrange': True},
                 hovermode='closest'
             )}, config)



@ app.callback(
    [Output('bubble_plot', 'figure'),
     Output('bubble_plot', 'config')],
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value'),
     ])
def update_bubble(sel_dataset, sel_standardized):
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    config = {'modeBarButtonsToRemove': plotly_buttons_remove, 'toImageButtonOptions': {
        'height': 500,
        'width': 700, 'scale': 5, 'filename': f'{dataset_name}_bubble.png'
    }}
    if sel_standardized == 'Original':
        if dict_plot_si[dataset_name] == 'No':
            return (empty_plot('Not enough datapoints' + '\n' + 'for a meaningful plot'), config)
    x_coords_list, y_coords_list, color_list, rv_id_list, scatter_size_list = unknown_essential_xy(
        dff)

    # make list of lines to add
    lines_add = []
    for x in [1.5, 2.5, 3.5, 4.5]:
        lines_add.append({'type': 'line',
                          'xref': 'x',
                          'yref': 'y',
                          'x0': x,
                          'y0': 0.5,
                          'x1': x,
                          'y1': 3.5,
                          'line': {'dash': 'dot',
                                   'color': 'grey'}
                          })
    for y in [1.5, 2.5]:
        lines_add.append({'type': 'line',
                          'xref': 'x',
                          'yref': 'y',
                          'x0': 0.5,
                          'y0': y,
                          'x1': 5.5,
                          'y1': y,
                          'line': {'dash': 'dot',
                                   'color': 'grey'}
                          })
    # print(scatter_size_list)
    return ({
        'data': [
            go.Scatter(
                x=x_coords_list,
                y=y_coords_list,
                text=rv_id_list,
                mode='markers',
                # opacity=1,
                hoverinfo='text',
                marker_size=scatter_size_list,
                marker={
                    # 'size': 15,
                    'line': {'width': 1.5, 'color': 'black'},
                    'color': color_list
                }
            )
        ],
        'layout': go.Layout(
            autosize=True,
            shapes=lines_add,
            # width=800,
            # height=500,
            xaxis=go.layout.XAxis(
                tickmode='array',
                tickvals=[1, 2, 3, 4, 5],
                ticktext=['least well<br>characterized', '', '', '',
                          'most well<br>characterized'],
                tickfont=dict(size=14),
                title='Annotation',
                showgrid=False

            ),
            yaxis=go.layout.YAxis(
                tickmode='array',
                tickvals=[1, 2, 3],
                ticktext=['non-essential', 'q-val < 0.05', 'q-val < 0.01'],
                tickangle=270,
                tickfont=dict(size=14),
                title='Essentiality',
                showgrid=False
            ),
            margin={'l': 30, 'b': 100, 't': 10, 'r': 10},
            legend={'x': 0, 'y': 1},
            hovermode='closest'
        )
    }, config)


@ app.callback(
    Output('dataset_metadata', 'children'),
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value')])
def print_dataset_metadata(sel_dataset, sel_standardized):
    if sel_standardized == 'Standardized':
        dataset_name = dict_si_to_std.get(sel_dataset, sel_dataset)
        dff = metadata[metadata['column_ID_std'] == dataset_name]
    else:
        dataset_name = dict_std_to_si.get(sel_dataset, sel_dataset)
        dff = metadata[metadata['column_ID_SI'] == dataset_name]

    # To handle nan data:
    dff = dff.fillna("No information available")

    text = [
        html.Strong('Summary'),
        html.Span(': ' + str(dff['meaning'].values[0])),
        html.Br(),
        html.Br(),
        html.Strong('Original Publication:'),
        html.Br(),
        html.A(dff['paper_title'].values[0],
               href=dff['paper_URL'].values[0], target='_blank'),
        html.Br(),
        html.Br(),
        html.Strong('Mtb strain'),
        html.Span(': ' + str(dff['Mtb strain'].values[0])),
        html.Br(),
        html.Br(),
        html.Strong('No of control replicates'),
        html.Span(': ' + str(dff['num replicates control'].values[0])),
        html.Br(),
        html.Br(),
        html.Strong('No of experimental replicates'),
        html.Span(': ' + str(dff['num replicates experimental'].values[0]))
    ]
    return text


@ app.callback(
    Output('sel_dataset_table', 'data'),
    [Input('sel_dataset', 'value'),
     Input('sel_standardized', 'value')])
def update_dataset_table(sel_dataset, sel_standardized):
    dff, _ = filter_dataset(sel_dataset, sel_standardized)
    dff = dff[['Rv_ID', 'gene_name', 'log2FC', 'q-val', 'id']]
    dff['q-val'] = dff['q-val'].astype(float).round(2)
    dff['log2FC'] = dff['log2FC'].astype(float).round(2)
    dff = dff.sort_values(by='log2FC')
    return dff.to_dict('records')


@ app.callback([Output('cog', 'figure'), Output('cog', 'config')],
               [Input('sel_dataset', 'value'),
                Input('sel_standardized', 'value'),
                Input('sel_cog', 'value'),
                Input('log2FC', 'value'),
                Input('q-val', 'value')])
def update_cog(sel_dataset, sel_standardized, sel_cog, log2FC, qval):
    dff, dataset_name = filter_dataset(sel_dataset, sel_standardized)
    config = {'modeBarButtonsToRemove': plotly_buttons_remove, 'toImageButtonOptions': {
        'height': 500,
        'width': 700, 'scale': 5, 'filename': f'{dataset_name}_bubble.png'
    }}
    if sel_standardized == 'Original':
        if dict_plot_si[dataset_name] == 'No':
            return (empty_plot('Not enough datapoints' + '\n' + 'for a meaningful plot'), config)
    if sel_cog == 'Under-represented':
        sel_subset_filter = (
            dff['q-val'] <= qval) & (dff['log2FC'] <= -log2FC)
        colorscale = 'Cividis'
    else:
        sel_subset_filter = (
            dff['q-val'] <= qval) & (dff['log2FC'] >= -log2FC)
        colorscale = 'Viridis'
    sel_subset = dff[sel_subset_filter]
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
    return ({'data': bar_data,
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
             }, config)


@ app.callback(
    Output('sel_gene_table', 'data'),
    [Input('sel_gene', 'value'),
     Input('sel_standardized_gene_table', 'value')])
def update_genes_table(selected_gene, sel_standardized_gene_table):
    if sel_standardized_gene_table == 'Standardized':
        dff = std_data.copy()
        metadata_col = 'column_ID_std'
    else:
        dff = si_data.copy()
        metadata_col = 'column_ID_SI'
    #if selected_gene in unique_Rvs:
    #    dff = dff[dff['Rv_ID'] == selected_gene]
    #elif selected_gene in unique_genes:
    #    dff = dff[dff['gene_name'] == selected_gene]
    if selected_gene in unique_Rvs_genes:
        name = selected_gene.split("/")[0]
        dff = dff[dff['Rv_ID'] == name]
    else:
        raise PreventUpdate
    metadata['paper'] = '[' + metadata['paper_title'] + \
        '](' + metadata['paper_URL'] + ')'
    metadata_cols_display = [col for col in metadata.columns if col not in [
        'plot_SI_graph', 'column_ID_std', 'column_ID_SI', 'paper_title', 'paper_URL']]
    metadata_trunc = metadata[[metadata_col] + metadata_cols_display]
    metadata_trunc = metadata_trunc.rename(columns={metadata_col: 'Expt'})
    merged_data = dff.merge(metadata_trunc, how='left', on='Expt')
    merged_data['q-val'] = np.round(merged_data['q-val'], 2)
    merged_data['log2FC'] = np.round(merged_data['log2FC'], 2)
    merged_data = merged_data.sort_values(by='q-val')
    if sel_standardized_gene_table == 'Standardized':
        merged_data['Expt'] = merged_data['Expt'].apply(
            lambda x: split_expt_name(x))
    return merged_data.to_dict('records')


@ app.callback(
    Output('gene_metadata', 'children'),
    [Input('sel_gene', 'value')])
def print_gene_metadata(sel_gene):
    #if sel_gene in unique_Rvs:
    #    sel_details = gene_metadata_df[gene_metadata_df['Rv_ID'] == sel_gene]
    #elif sel_gene in unique_genes:
    #    sel_details = gene_metadata_df[gene_metadata_df['gene_name'] == sel_gene]
    #    sel_details = si_data[si_data['gene_name'] == sel_gene]
    if sel_gene in unique_Rvs_genes:
        sel_details = gene_metadata_df[gene_metadata_df['normalized_name'] == sel_gene]
    else:
        return "gene not found"
    text = [
        html.Span(list(sel_details['Description'])[0]),
        html.Br(),
        html.Strong('mBio Call: '),
        html.Span(list(sel_details['Final Call'])[0]),
        html.Br(),
        html.Strong('Tuberculist functional category: '),
        html.Span(list(sel_details['tuberculist_category'])[0])
    ]
    return text

# Co-essentiality

@ app.callback(
    Output('correlation_plot', 'figure'),  # Corrected the ID to match the dcc.Graph component
    [Input('sel_gene', 'value'),
    Input('sel_warped_gene', 'value')]
)

def correlation_plot_query(sel_gene, sel_warped_gene):

    sel_gene = sel_gene.split("/")[0]
    
    list_rvid_NN1, list_rvid_NN2 = get_NN12(sel_gene, df_interact)
    
    list_rvid_x = list_rvid_NN2.copy()
    list_rvid_y = list_rvid_NN2.copy()
    
    # Create a set of significant interactions for quick lookup
    significant_pairs = set(df_interact[['lead_gene', 'partner_gene']].itertuples(index=False, name=None))
    
    # Create a dictionary to map gene pairs to p-values
    p_value_dict = {(row['lead_gene'], row['partner_gene']): row['p_value_FDR'] 
                    for index, row in df_interact.iterrows()}


    # Cut the total plots to 10 in both axis to ensure the plot is shown                
    if len(list_rvid_x) > 10 and len(list_rvid_y) > 10:
        common_sample = sample(list_rvid_x[1:], 9)
        list_rvid_x = [sel_gene] + common_sample
        list_rvid_y = [sel_gene] + common_sample
        #list_rvid_y = [sel_gene] + sample(list_rvid_y[1:], 9)
        title = f"Showing 10 Randomized genes that have a significant co-essentiality correlation with {sel_gene} (first and second degree correlations)"
    else: 
        title = f"Showing all genes that have a significant co-essentiality correlation with {sel_gene} (first and second degree correlations) "

    list_gene_names_x = [f"{dict_rvid_to_name[rvid]}<br>({rvid})"  if dict_rvid_to_name[rvid] != rvid else f"<br>{rvid}" for rvid in list_rvid_x]
    list_gene_names_y = [f"{dict_rvid_to_name[rvid]}<br>({rvid})"  if dict_rvid_to_name[rvid] != rvid else f"<br>{rvid}" for rvid in list_rvid_y]
   
    # Making the subplots
     
    fig = make_subplots(rows=len(list_rvid_x), cols=len(list_rvid_y)) #subplot_titles=[f'{var1} vs {var2}' for var1 in list_rvid for var2 in list_rvid])        
    
    for i in range(len(list_rvid_x)):
        for j in range(len(list_rvid_y)):
            x_rvid = list_rvid_x[i]
            y_rvid = list_rvid_y[j]

            if sel_warped_gene == "Warped":
                x = get_warped_screen_for_gene(x_rvid, df_lfc, warped_screens) # maybe is better to have a dataframe with this values instead of calculating on the fly
                y = get_warped_screen_for_gene(y_rvid, df_lfc, warped_screens)
            else:
                x = df_lfc.loc[list_rvid_x[i]].values
                y = df_lfc.loc[list_rvid_y[j]].values
            # Check if the current pair is in the significant interactions
            is_significant = (x_rvid, y_rvid) in significant_pairs or (y_rvid, x_rvid) in significant_pairs
            # Set alpha value depending on significance
            alpha_val = 0.9 if is_significant else 0.2
            
            
            # Creating internal reference of the subplots an
            xref_int = ((i*len(list_rvid_y))+j)+1
            yref_int = ((i*len(list_rvid_y))+j)+1

            fig.add_trace(go.Scatter(x=x, y=y, mode='markers',opacity=alpha_val, showlegend=False), row=i+1, col=j+1) # Plus one because plottly uses the starting 1 convention ### row=i+1, col=j+1
            
            
            if is_significant:
                # Retrieve the p-value for the current gene pair using both orientations
                # Using the `get` method to avoid KeyError if the pair is not found
                p_value = p_value_dict.get((list_rvid_x[i], list_rvid_y[j]), 
                                           p_value_dict.get((list_rvid_y[j], list_rvid_x[i]), None))

                # Format the p-value in scientific notation
                p_value_text = f'p={p_value:.2e}'
                
                # Assuming a standard subplot size, calculate annotation position
                fig.add_annotation(xref=f'x{xref_int}', 
                                   yref=f'y{yref_int}',
                                   x=0.5, y=y.max()+1, showarrow=False,
                                   text=f'p={p_value:.2e}', font=dict(size=10),
                                   xanchor='center', yanchor='bottom')

# Making the Annotations of the plots with the p_values
    for i, rvid in enumerate(list_gene_names_x):
        fig.update_yaxes(title=list_gene_names_x[i], row=i+1, col=1)

    # Update x-axis labels for the bottom row only
    for j, rvid in enumerate(list_gene_names_y):
        fig.update_xaxes(title=list_gene_names_y[j], row=len(list_rvid_y), col=j+1)
        

    fig.update_layout(title=title, width=1400, height=1400, hovermode='closest')

    return fig

@ app.callback(
    Output('sel_coesen_table', 'data'),
    [Input('sel_gene', 'value')])

def update_coesen_table(sel_gene):
    sel_gene = sel_gene.split("/")[0]
    list_rvid_NN1, list_rvid_NN2 = get_NN12(sel_gene, df_interact)
    #dff = df_interact[((df_interact["lead_gene"] == sel_gene) | (df_interact["partner_gene"] == sel_gene))]
    dff = df_interact[ (df_interact.lead_gene.isin(list_rvid_NN1)) | (df_interact.partner_gene.isin(list_rvid_NN1))].copy()
    dff["p_value_FDR"] = dff["p_value_FDR"].apply(lambda x: f"{x:.2e}")
    dff['lead_gene'] = dff['lead_gene'].map(lambda x: f"{dict_rvid_to_name.get(x, x)} ({x})" if dict_rvid_to_name.get(x, x) != x else x)
    dff['partner_gene'] = dff['partner_gene'].map(lambda x: f"{dict_rvid_to_name.get(x, x)} ({x})" if dict_rvid_to_name.get(x, x) != x else x)
    return dff.to_dict('records')


if __name__ == '__main__':
    app.run_server(debug=True)
