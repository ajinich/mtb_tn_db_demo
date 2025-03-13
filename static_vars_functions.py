import os
import urllib.parse
from io import StringIO
from random import random

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
import json
from explorer_lib import get_NN12

external_stylesheets = [dbc.themes.UNITED]
# path_data = '../../data/'

# read in data
#std_data = pd.read_csv('https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/standardized_data_dash.tsv',
#                       sep='\t', dtype={'Rv_ID': str, 'gene_name': str, 'Description': str, 'Expt': str, 'log2FC': float, 'q-val': float})

std_data = pd.read_csv('data/standardized_data_dash.tsv',
                       sep='\t', dtype={'Rv_ID': str, 'gene_name': str, 'Description': str, 'Expt': str, 'log2FC': float, 'q-val': float})
si_data = pd.read_csv('data/si_data_dash.tsv', sep='\t',
                      dtype={'Rv_ID': str, 'gene_name': str, 'Description': str, 'Expt': str, 'log2FC': float, 'q-val': float})
#si_data = pd.read_csv('https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/si_data_dash.tsv', sep='\t',
#                      dtype={'Rv_ID': str, 'gene_name': str, 'Description': str, 'Expt': str, 'log2FC': float, 'q-val': float})

#gene_metadata_df = pd.read_csv(
#    'https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/gene_metadata_dash.tsv', sep='\t')

gene_metadata_df = pd.read_csv(
    'data/gene_metadata_dash.tsv', sep='\t')

#metadata = pd.read_csv(
#    'https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/col_desc_dash.tsv', sep='\t')

metadata = pd.read_csv(
    'data/col_desc_dash.tsv', sep='\t')

# make dictionaries that lookup si_data from std_data and vice versa
dict_std_to_si = dict(zip(metadata.column_ID_std, metadata.column_ID_SI))
dict_si_to_std = dict(zip(metadata.column_ID_SI, metadata.column_ID_std))

# make dictionary that tells you if SI data has enough datapoints to be plotted
dict_plot_si = dict(zip(metadata.column_ID_SI, metadata.plot_SI_graph))

# make lists of unique expts/rvids/genes for dropdown lists
unique_expts = list(metadata['column_ID_std'].unique()) + \
    list(metadata['column_ID_SI'].unique())
unique_expts = [x for x in unique_expts if str(x) != 'nan']
unique_filter_col = ["All", "in vitro/cell/in vivo", "in vitro media", "carbon source", "stress", "Microarray or TnSeq", "mouse strain", "cell type", "Mtb strain"]


#unique_Rvs = sorted(gene_metadata_df.Rv_ID.to_list())
#unique_genes = sorted(gene_metadata_df.gene_name.to_list())
#unique_genes = [x for x in unique_genes if x != '-']
unique_Rvs_genes =  sorted(gene_metadata_df["normalized_name"].to_list())

co_genes = [x for x in unique_Rvs_genes if (gene_metadata_df[gene_metadata_df["Rv_ID"]==x.split("/")[0]]["coessential_signal"] == "Positive").any()]

# do data wrangling on si and std data for dash requirements
std_data['id'] = std_data['Rv_ID']
std_data.set_index('id', inplace=True, drop=False)
si_data['id'] = si_data['Rv_ID']
si_data.set_index('id', inplace=True, drop=False)

# read in cog annotations, cogs_desc should be read as a pd.Series
# path_annotation = '../../data/annotations/'
#cogs_df = pd.read_csv(
#    'https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/annotations/all_cogs.csv')

cogs_df = pd.read_csv(
    'data/annotations/all_cogs.csv')
#cogs_desc = pd.read_csv(
#    'https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/annotations/cog_names.csv', header=0, index_col=0, squeeze=True)
cogs_desc = pd.read_csv(
    'data/annotations/cog_names.csv', header=0, index_col=0, squeeze=True)


# select columns of dataset metadata to show in analyse genes table
metadata_cols_display = [col for col in metadata.columns if col not in [
    'plot_SI_graph', 'column_ID_std', 'column_ID_SI', 'paper_title', 'paper_URL']]
sel_gene_table_columns = [{"name": i,
                           "id": i,
                           } for i in [
    'Rv_ID', 'gene_name', 'Expt', 'log2FC', 'q-val'] + metadata_cols_display]

# show paper column as markdown so that html links work
sel_gene_table_columns.append(
    {"name": 'paper', "id": 'paper', "presentation": 'markdown'})

# show paper column as markdown so that html links work for track view
sel_gene_table_columns.append(
    {"name": 'Track     View', "id": 'Track View', "presentation": 'markdown'})

# list plotly buttons to remove from all graphs
plotly_buttons_remove = [
    'pan2d', 'lasso2d', 'select2d', 'hoverClosestCartesian', 'hoverCompareCartesian',
    'zoomIn2d', 'zoomOut2d', 'toggleSpikelines', 'autoScale2d', 'zoom2d'
]

# Adding information for co-essentiality

#df_interact = pd.read_csv("https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/coessential/interaction_pairs.csv")
#df_lfc = pd.read_csv("https://raw.githubusercontent.com/ajinich/mtb_tn_db_demo/master/data/coessential/values.csv", index_col="Rv_ID")

df_interact = pd.read_csv("data/coessential/interaction_pairs.csv")
df_lfc = pd.read_csv("data/coessential/values.csv", index_col="Rv_ID")

with open('data/coessential/rvid_to_name.json', 'r') as json_file:
    dict_rvid_to_name = json.load(json_file)

sel_coesen_table_columns = columns = [
    {"name": "Lead Gene", "id": "lead_gene"},
    {"name": "Partner Gene", "id": "partner_gene"},
    {"name": "P-Value FDR", "id": "p_value_FDR"},
    {"name": "Sorted Gene Pair", "id": "sorted_gene_pair"},
    {"name": "Outlier Driven Flag", "id": "outlier_driven_flag"}
]

data_warped_screens = np.load('data/coessential/warped_screens.npz')
warped_screens = tuple(data_warped_screens[key] for key in data_warped_screens)

# Adding data for trackviews
track_df = pd.read_csv("data/track_view_urls.tsv.gz", sep="\t", compression="gzip")


def empty_plot(label_annotation):
    """Returns a dictionary with elements for an empty plot with centered text

    Args:
        label_annotation (str): Text to display

    Returns:
        dict: Elements for plotting
    """
    return {
        "layout": {
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "annotations": [
                {"text": label_annotation,
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {"size": 16}
                 }]
        }}


def split_expt_name(expt_name):
    expt_name_split = expt_name.split('_vs_')
    new_expt_name = expt_name_split[0] + '_vs_' + '\n' + expt_name_split[1]
    return new_expt_name


def unknown_essential_xy(selected_data):
    """Generate plotting data for bubble plot

    Args:
        selected_data (pd.DataFrame): Data filtered on selected experiment

    Returns:
        Multiple outputs required for plotting
    """

    df_vis = selected_data.dropna(subset=['annotation_score'])
    df_vis = df_vis.reset_index(drop=True)

    # NO MORE RANDOMIZATION! Instead this:
    rv_id_list = []
    x_coords_list = []
    y_coords_list = []
    color_list = []
    scatter_size_list = []

    for ann in [1, 2, 3, 4, 5]:
        for qq in [1, 2, 3]:

            df = df_vis[(df_vis.q_val_D.values == qq) &
                        (df_vis.annotation_score.values == ann)]
            # update rv_id list
            rv_id_list += list(df.Rv_ID.values)

            if df.shape[0] < 30:
                scatter_size = 6
                edge_param = 0.40
            elif df.shape[0] < 100:
                scatter_size = 6
                edge_param = 0.40
            elif df.shape[0] < 200:
                scatter_size = 3
                edge_param = 0.45
            elif df.shape[0] < 400:
                scatter_size = 3
                edge_param = 0.45
            else:
                scatter_size = 3
                edge_param = 0.45

            # Update scatter marker size
            scatter_size_list += [scatter_size] * df.shape[0]

            if ann <= 2 and qq >= 2:
                color_temp = '#2b7bba'
            else:
                color_temp = '#585858'

            # update color_list
            color_list += [color_temp] * df.shape[0]

            # num_sqrs = int(np.ceil(np.sqrt(df.shape[0])))
            num_sqrs = np.max([int(np.ceil(np.sqrt(df.shape[0]))), 10])

            xrange = np.linspace(ann - edge_param, ann + edge_param, num_sqrs)
            yrange = np.linspace(qq + edge_param, qq - edge_param, num_sqrs)

            coords = list(itertools.product(yrange, xrange))
            coords = coords[: df.shape[0]]
            x_coords = [c[1] for c in coords]
            y_coords = [c[0] for c in coords]
            # update coordinates
            x_coords_list += x_coords
            y_coords_list += y_coords

    return x_coords_list, y_coords_list, color_list, rv_id_list, scatter_size_list


def filter_dataset(sel_dataset, sel_standardized):
    """Using outputs from sel_dataset and sel_standardized return appropriate data and name

    Args:
        sel_dataset (str): Selected dataset eg: griffin_glycerol_vs_mbio_H37Rv
        sel_standardized (str): Selected standardized eg: 'Original'

    Returns:
        pd.DataFrame: Filtered data
        str: Dataset name.
    """
    if sel_standardized == 'Standardized':
        dataset_name = dict_si_to_std.get(sel_dataset, sel_dataset)
        dff = std_data[std_data['Expt'] == dataset_name]
    else:
        dataset_name = dict_std_to_si.get(sel_dataset, sel_dataset)
        dff = si_data[si_data['Expt'] == dataset_name]
    return dff, dataset_name


def update_download_dataset(dff, dataset_name):
    """Given data filtered by experiment, provide data for downloading

    Args:
        dff (pd.DataFrame): Data filtered by experiment
        dataset_name (str): Dataset name

    Returns:
        str: download string
        str: download file name
    """
    dff = dff.copy(deep=True)
    dff.reset_index(inplace=True, drop=True)
    dff = dff.drop(columns='id')
    csv_string = dff.to_csv(encoding='utf-8', sep='\t', index=False)
    csv_string = "data:text/plain;charset=utf-8," + \
        urllib.parse.quote(csv_string)
    download_string = dataset_name + '.tsv'
    return csv_string, download_string


# TODO: decide what to do with NAs
def update_num_significant(dff, log2FC, qval):
    """Given filtered data, log2FC and qval return text indicating num significant genes

    Args:
        dff (pd.DataFrame): Data filtered by expt
        log2FC (float): log2FC cutoff
        qval (float): qval cutoff

    Returns:
        list: List of html elements as text
    """
    dff = dff.dropna(axis=0, subset=['log2FC', 'q-val'])
    num_neg_sig = dff[(dff['q-val'] <= qval) &
                      (dff['log2FC'] <= -log2FC)].shape[0]
    num_pos_sig = dff[(dff['q-val'] <= qval) &
                      (dff['log2FC'] >= log2FC)].shape[0]
    text = [html.Span(f'Number of neg sig : {num_neg_sig} ; Number of pos sig : {num_pos_sig}'),
            ]
    return text

def get_warped_screen_for_gene(gene_id, df, warped_screens):
    """
    Retrieve the warped screen data for a specific gene.

    Parameters:
    gene_id (str): The identifier of the gene.
    df (DataFrame): The DataFrame containing gene information.
    warped_screens (ndarray): The array of warped screens data.

    Returns:
    ndarray: The warped screen data corresponding to the given gene_id.
    """
    # Get the index of the gene in the DataFrame
    gene_index = df.index.get_loc(gene_id)
    
    # Extract the corresponding row from warped_screens
    gene_warped_screen = warped_screens[gene_index]

    return gene_warped_screen