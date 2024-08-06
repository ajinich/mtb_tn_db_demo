import pandas as pd

import os
import re

import matplotlib.pyplot as plt 
import numpy as np

# from sklearn.preprocessing import normalize
# from sklearn.metrics import pairwise
# import torch

# import holoviews as hv 
# import bokeh.io
# import bokeh.plotting
# import colorcet as cc

def import_interact_lfc_uniprot(norm, fn_interact, fn_lfc, fn_annot):
    """
    Imports and preprocesses interaction and log-fold change (lfc) data, along with gene annotations.

    Parameters:
    norm (bool): If True, normalizes the lfc data using L2 normalization. If False, uses raw data.
    fn_interact (str): Filename or path to the Excel file containing gene interaction data.
    fn_lfc (str): Filename or path to the CSV file containing log-fold change data for genes.
    fn_annot (str): Filename or path to the Excel file containing gene annotations.

    Returns:
    tuple: A tuple containing four elements:
        - df_interact (DataFrame): DataFrame of gene interaction data.
        - df_lfc (DataFrame): DataFrame of normalized (if norm is True) or raw lfc data.
        - df_annot (DataFrame): DataFrame of gene annotations.
        - dict_rvid_to_name (dict): Dictionary mapping Rv IDs to gene names.

    The function performs the following operations:
    1. Reads gene interaction data from an Excel file.
    2. Reads and optionally normalizes lfc data from a CSV file.
    3. Reads gene annotations from an Excel file, processes Rv IDs, and gene names.
    4. Creates a dictionary for mapping Rv IDs to gene names for easy lookup.

    Note:
    - The lfc data CSV is expected to have genes as rows and experimental conditions as columns.
    - Gene annotations should include 'Gene names' from which Rv IDs and gene names are extracted.
    """
    
    df_interact = pd.read_excel(fn_interact)

    df_lfc = pd.read_csv(fn_lfc, index_col=0)
    df_lfc.dropna(axis=0, inplace=True)

    cols_data = df_lfc.columns[1:]

    if norm: 
        X = df_lfc[cols_data].values
        X_norm = normalize(X, norm='l2', axis=0)
        df_lfc_norm = df_lfc.copy()
        df_lfc_norm[cols_data] = X_norm
        df_lfc = df_lfc_norm.copy()
        
    else:
        df_lfc = df_lfc.copy()

    # loading annotations: 
    df_annot = pd.read_excel(fn_annot)
    df_annot = df_annot.fillna('')

    re_str = 'Rv\d\d\d\dc?'
    list_rvids = [re.findall(re_str, str_temp)[0] for str_temp in df_annot['Gene names']]
    df_annot['Rv_ID'] = list_rvids

    list_gene_names = [gn.split()[0] for gn in df_annot["Gene names"]]
    df_annot['gene_names'] = list_gene_names

    df_rvid_to_name = df_annot[['Rv_ID', 'gene_names']].copy() 

    dict_rvid_to_name = {}
    for index, row in df_rvid_to_name.iterrows():
        dict_rvid_to_name[row.Rv_ID] = row.gene_names


    return df_interact, df_lfc, df_annot, dict_rvid_to_name


def get_NN12(rvid_query, df_interact):
    """
    Generate lists of first and second nearest neighbors (NN1 and NN2) for a given gene.

    Parameters:
    rvid_query (str): The RvID of the seed gene for which nearest neighbors are to be found.
    df_interact (DataFrame): A DataFrame containing gene interaction pairs (columns 'lead_gene' and 'partner_gene').

    Returns:
    tuple: A tuple containing two lists:
        - The first list contains the first nearest neighbors (NN1) of the query gene.
        - The second list contains the second nearest neighbors (NN2) of the query gene.

    The function works as follows:
    1. NN1 Identification:
       - Finds all genes that are directly correlated with the query gene (either as a lead_gene or partner_gene).
       - These genes are considered the first nearest neighbors (NN1).

    2. NN2 Identification:
       - Extends the search to genes that are correlated with the genes in the NN1 list.
       - These additional genes are considered the second nearest neighbors (NN2).

    Both lists are deduplicated and sorted to facilitate further analysis.
    """
    # first nearest neighbors: 
    df_NN1 = df_interact[(df_interact.lead_gene==rvid_query) | (df_interact.partner_gene==rvid_query)].copy()
    list_rvid_NN1 = list(set(df_NN1.lead_gene.tolist() + df_NN1.partner_gene.tolist()))
    list_rvid_NN1.sort()

    # second nearest neighbors: 
    df_NN2 = df_interact[ (df_interact.lead_gene.isin(list_rvid_NN1)) | (df_interact.partner_gene.isin(list_rvid_NN1))].copy()
    list_rvid_NN2 = list(set(df_NN2.lead_gene.tolist() + df_NN2.partner_gene.tolist()))
    list_rvid_NN2.sort()
    
    return list_rvid_NN1, list_rvid_NN2


def correlation_tile_plot(df_lfc, df_interact, list_rvid_x, list_rvid_y, fig_size, cols, dict_rvid_to_name, list_subset=[], gene_names = True ):
    
    if gene_names:
        list_gene_names_x = [dict_rvid_to_name[rvid] if rvid in dict_rvid_to_name.keys() else rvid for rvid in list_rvid_x]
        list_gene_names_y = [dict_rvid_to_name[rvid] if rvid in dict_rvid_to_name.keys() else rvid for rvid in list_rvid_y]
    else:
        list_gene_names_x = list_rvid_x
        list_gene_names_y = list_rvid_y
    
    # Create a set of significant interactions for quick lookup
    significant_pairs = set(df_interact[['lead_gene', 'partner_gene']].itertuples(index=False, name=None))
    
    # Create a dictionary to map gene pairs to p-values
    p_value_dict = {(row['lead_gene'], row['partner_gene']): row['p_value_FDR'] 
                    for index, row in df_interact.iterrows()}
    
    fig, axs = plt.subplots(len(list_rvid_x), len(list_rvid_y), figsize=fig_size)
    if max([len(list_rvid_x), len(list_rvid_y)]) >= 10:
        FontSize = 20
        size_param = 40
    elif max([len(list_rvid_x), len(list_rvid_y)]) <= 5:
        FontSize = 20
        size_param = 60
    else:
        FontSize = 20
        size_param = 60
        
    
    for i in range(len(list_rvid_x)):
        for j in range(len(list_rvid_y)):
            x_rvid = list_rvid_x[i]
            y_rvid = list_rvid_y[j] 

            x = df_lfc[df_lfc.Rv_ID==x_rvid].values[0][1:]
            y = df_lfc[df_lfc.Rv_ID==y_rvid].values[0][1:]
            
            # Check if the current pair is in the significant interactions
            is_significant = (x_rvid, y_rvid) in significant_pairs or (y_rvid, x_rvid) in significant_pairs
            # Set alpha value depending on significance
            alpha_val = 0.7 if is_significant else 0.05

            axs[i,j].scatter(x, y, s = size_param, alpha = alpha_val, edgecolors='k', linewidths=1)
            
            if is_significant:
                # Retrieve the p-value for the current gene pair using both orientations
                # Using the `get` method to avoid KeyError if the pair is not found
                p_value = p_value_dict.get((list_rvid_x[i], list_rvid_y[j]), 
                                           p_value_dict.get((list_rvid_y[j], list_rvid_x[i]), None))

                # Format the p-value in scientific notation
                p_value_text = f'p={p_value:.2e}'
                # Annotate the scatter plot with the p-value
#                 axs[i, j].text(0.2, 0.9, p_value_text, fontsize=9, ha='center', va='center', transform=axs[i, j].transAxes)
                axs[i,j].set_title(p_value_text, fontsize = 10)

            if i == len(list_rvid_x)-1:
                axs[i,j].set_xlabel(list_gene_names_x[j], fontsize = FontSize)
            
            if j == 0:
                axs[i,j].set_ylabel(list_gene_names_y[i], fontsize = FontSize)
            
            if list_rvid_x == list_rvid_y:
                if x_rvid in list_subset or y_rvid in list_subset:
                    axs[i,j].set_facecolor('xkcd:lightblue')
                if x_rvid in list_subset and y_rvid in list_subset:
                    axs[i,j].set_facecolor(cols[-2])
                if i==j:
                    axs[i,j].set_facecolor(cols[-3])


def interactive_scatter_grid(df_lfc, list_rvid):

    cols_data = df_lfc.columns[1:].tolist()
    
    df_xy = df_lfc[ df_lfc.Rv_ID.isin(list_rvid) ][ ['Rv_ID']+cols_data ].copy()
    df_xy = df_xy.set_index('Rv_ID').T.rename_axis('screen').reset_index()

    list_hv = []
    for i in range(len(list_rvid)):
        hv_temp = df_xy.hvplot.scatter(x = list_rvid[i], y = list_rvid, width = 300, height = 300, size = 200, line_color='k', 
                                       line_width=3, hover_cols = ['screen'], subplots=True, fontsize = {'xlabel': '15pt'}, xlabel = list_rvid[i] ).cols(len(list_rvid))
        list_hv.append(hv_temp)

    return list_hv



def load_ESM_embeddings_and_df():
    # path to Mtb/M.smeg proteome: 
    fn = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/df_mtb_smeg_umap.csv'
    df_mtb_smeg = pd.read_csv(fn)
    df_mtb_smeg = df_mtb_smeg.fillna('')

    list_rvids = []
    re_str = 'Rv\d\d\d\dc?'
    for str_temp in df_mtb_smeg['Gene names']:
        re_match = re.findall(re_str, str_temp)
        if len(re_match):
            list_rvids.append(re_match[0])
        else:
            list_rvids.append('')
    df_mtb_smeg['Rv_ID'] = list_rvids

    # Path to ESM embeddings / representations: 
    path_rep = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/up_mtb_smeg_reprs/'
    EMB_LAYER = 33

    Xs_list = []
    list_err = []
    list_entries = df_mtb_smeg.Entry.tolist()

    for entry in list_entries:
        fn_full = os.path.join(path_rep, entry+'.pt')
        try:
            embs = torch.load(fn_full)
            Xs_list.append(embs['mean_representations'][EMB_LAYER])
        except:
            list_err.append(entry)
    X = torch.stack(Xs_list, dim=0).numpy()

    return X, df_mtb_smeg, list_entries 


def load_ESM_embeddings_and_df_10_prots():

    path_df = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/sakila_ESM/Proteomes/'
    list_fn = [os.path.join(path_df, fn) for fn in os.listdir(path_df)]
    list_df = [pd.read_csv(fn, sep='\t') for fn in list_fn]
    df_orgs = pd.concat(list_df, axis = 0)

    col = 'Gene names'
    df_orgs[col] = df_orgs[col].fillna('')
    list_rvids = []
    re_str = 'Rv\d\d\d\dc?'
    for str_temp in df_orgs[col]:
        re_match = re.findall(re_str, str_temp)
        if len(re_match):
            list_rvids.append(re_match[0])
        else:
            list_rvids.append('')
    df_orgs['Rv_ID'] = list_rvids

    path_rep = '/home/ajinich/Dropbox/KyuRhee/unknown_function/unknown_redox/data/mohammed/sakila_ESM/TotalProteomes/'
    EMB_LAYER = 33

    list_entries = [fn.split('.')[0] for fn in os.listdir(path_rep)]
    df_orgs = df_orgs[df_orgs.Entry.isin(list_entries)]
    df_orgs.reset_index(inplace=True, drop = True)

    Xs_list = []
    list_err = []
    list_entries = df_orgs.Entry.tolist()

    for entry in list_entries:
        fn_full = os.path.join(path_rep, entry+'.pt')
        try:
            embs = torch.load(fn_full)
            Xs_list.append(embs['mean_representations'][EMB_LAYER])
        except:
            list_err.append(entry)
    X = torch.stack(Xs_list, dim=0).numpy()

    return X, df_orgs, list_entries


def get_similar_prots(X, df, ind_query, list_entries, perc_th):

    x_test = X[ind_query,:]
    x_test = x_test.reshape(1, -1)

    mat_sim = pairwise.cosine_similarity(x_test, X)
    mat_d = 1 - mat_sim
    d_th = np.percentile(mat_d, perc_th)
    ind_th = list(np.where(mat_d <= d_th)[1])

    list_entries_th = [list_entries[i] for i in ind_th]

    df_th = pd.DataFrame()
    df_th['Entry'] = list_entries_th
    df_th['d']= mat_d[0, ind_th]
    df_th.sort_values(by = 'd', inplace=True)

    df_th = df_th.merge(df, on = 'Entry', how = 'left')

    return df_th