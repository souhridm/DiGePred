#!/usr/bin/env python

import pandas as pd  # version 0.25.1
import numpy as np
import _pickle as pickle
import networkx as nx  # version 1.9
import argparse
import itertools
import datetime

now = datetime.datetime.now()
month = str(now.strftime("%m"))
day = str(now.strftime("%d"))
year = str(now.strftime("%Y"))
hour = str(now.strftime("%H"))
minute = str(now.strftime("%M"))

## Load pathway data files

reactome_gene_to_path_codes = pickle.load(open('~/DiGePred/data/pathways/path-to-{reactome_gene_to_path_codes}-file', 'rb'))
reactome_path_to_genes = pickle.load(open('~/DiGePred/data/pathways/path-to-{reactome_path_to_genes}-file', 'rb'))
reactome_path_code_to_name = pickle.load(open('~/DiGePred/data/pathways/path-to-{reactome_path_code_to_name}-file', 'rb'))
reactome_tot_codes_in_path = pickle.load(open('~/DiGePred/data/pathways/path-to-{reactome_tot_codes_in_path}-file', 'rb'))

kegg_list_genes = pickle.load(open('~/DiGePred/data/pathways/path-to-{kegg_list_genes}-file', 'rb'))
kegg_full_gene_name = pickle.load(open('~/DiGePred/data/pathways/path-to-{kegg_full_gene_name}-file', 'rb'))
kegg_gene_name = pickle.load(open('~/DiGePred/data/pathways/path-to-{kegg_gene_name}-file', 'rb'))
kegg_gene_to_path_codes = pickle.load(open('~/DiGePred/data/pathways/path-to-{kegg_gene_to_path_codes}-file', 'rb'))
kegg_path_to_genes = pickle.load(open('~/DiGePred/data/pathways/path-to-{kegg_path_to_genes}-file', 'rb'))
kegg_path_code_to_name = pickle.load(open('~/DiGePred/data/pathways/path-to-{kegg_path_code_to_name}-file', 'rb'))

## Load phenotype data files

hpo_full_name_to_codes = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{hpo_full_name_to_codes}-file', 'rb'))
hpo_code_to_all_names = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{hpo_code_to_all_names }-file', 'rb'))
similar_codes = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{similar_codes}-file', 'rb'))
hpo_gene_to_code = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{hpo_gene_to_code}-file', 'rb'))
hpo_code_to_gene = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{hpo_code_to_gene}-file', 'rb'))
hpo_name_to_code = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{hpo_name_to_code}-file', 'rb'))
hpo_code_to_name = pickle.load(open('~/DiGePred/data/phenotypes/path-to-{hpo_code_to_name}-file', 'rb'))

## Load co-expression data files

coexpress_dict = pickle.load(open('~/DiGePred/data/coex/path-to-{coexpress_dict}-file', 'rb'))

## Load network data files

G_ppi = nx.read_dot('~/DiGePred/data/networks/path-to-{UCSC_ppi_dot}-file')
G_pwy = nx.read_dot('~/DiGePred/data/networks/path-to-{UCSC_pwy_dot}-file')
G_txt = nx.read_dot('~/DiGePred/data/networks/path-to-{UCSC_txt-dot}-file')

dists_ppi = pickle.load(open('~/DiGePred/data/networks/path-to-{ppi_dists_dict}-file', 'rb'))
dists_pwy = pickle.load(open('~/DiGePred/data/networks/path-to-{pwy_dists_dict}-file', 'rb'))
dists_txt = pickle.load(open('~/DiGePred/data/networks/path-to-{txt_dists_dict}-file', 'rb'))

## Load evoltuonary biology and genomics feature data files

lof_dict = pickle.load(open('~/DiGePred/data/networks/path-to-{lof_pli_dict}-file', 'rb'))
hap_insuf_dict = pickle.load(open('~/DiGePred/data/networks/path-to-{happloinsufficiency_dict}-file', 'rb'))
protein_age_dict = pickle.load(open('~/DiGePred/data/networks/path-to-{protein_age_dict}-file', 'rb'))
dNdS_full_dict = pickle.load(open('~/DiGePred/data/networks/path-to-{dNdS_full_dict}-file', 'rb'))
dNdS_avg_dict = pickle.load(open('~/DiGePred/data/networks/path-to-{dNdS_avg_dict}-file', 'rb'))
gene_ess_dict = pickle.load(open('~/DiGePred/data/networks/path-to-{Gene_Essentiality_dict}-file', 'rb'))

parser = argparse.ArgumentParser(description='Get feature dataframe')

parser.add_argument('-g', '--genes', type=str,
                    help='genes to get feature dataframe',
                    dest='genes', required=False,
                    metavar='')

parser.add_argument('-p', '--pairs', type=str,
                    help='pairs to get feature dataframe',
                    dest='pairs', required=False,
                    metavar='')

parser.add_argument('-m', '--model', type=str, default='unaffected no gene overlap',
                    help='model',
                    dest='model', required=False,
                    metavar='')

parser.add_argument('-n', '--name', type=str, default='{y}-{m}-{d}-{hr}{mi}'.format(m=month, d=day, y=year, hr=hour, mi=minute),
                    help='name',
                    dest='name', required=False,
                    metavar='')

genes = []
pairs = []

if parser.parse_args().pairs:
    pairs_list = open(parser.parse_args().pairs).read().split('\n')[:-1]

    for line in pairs_list:
        g1 = line.split(',')[0].strip().rstrip()
        g2 = line.split(',')[1].strip().rstrip()

        pairs.append(tuple(sorted([g1, g2])))

    pairs = sorted(set(pairs))

elif parser.parse_args().genes:

    genes = sorted(set(open(parser.parse_args().genes).read().split('\n')[:-1]))
    pairs = itertools.combinations(set(genes), 2)

else:
    print('No input file found!')

project_name = parser.parse_args().name
model = parser.parse_args().model


## Load DiGePred models

clfs = dict()

clfs['permuted'] = pd.read_pickle('~/DiGePred/models/path-to-{permuted_model}-file')
clfs['random'] = pd.read_pickle('~/DiGePred/models/path-to-{random_model}-file')
clfs['matched'] = pd.read_pickle('~/DiGePred/models/path-to-{matched_model}-file')
clfs['unaffected'] = pd.read_pickle('~/DiGePred/models/path-to-{unaffected_model}-file')
clfs['all-digenic-vs-unaffected'] = pd.read_pickle('~/DiGePred/models/path-to-{all-digenic-vs-unaffected_model}-file')
clfs['unaffected-no-gene-overlap'] = pd.read_pickle('~/DiGePred/models/path-to-{unaffected-no-gene-overlap_model}-file')
clfs['random-no-gene-overlap'] = pd.read_pickle('~/DiGePred/models/path-to-{random-no-gene-overlap_model}-file')

## Function to get feature values as a pandas dataframe
def get_features(pairs):

    new_list_pairs = [tuple(sorted(p)) for p in list(set(pairs))]

    all_data = []

    for x in new_list_pairs:

        data = np.empty((1, 21))

        #  Pathway
        path1 = []
        path2 = []

        if x[0] in kegg_gene_to_path_codes or x[0] in reactome_gene_to_path_codes:

            if x[0] in kegg_gene_to_path_codes:
                path1 = kegg_gene_to_path_codes[x[0]]

            if x[0] in reactome_gene_to_path_codes:
                path1.extend(reactome_gene_to_path_codes[x[0]])

        if x[1] in kegg_gene_to_path_codes or x[1] in reactome_gene_to_path_codes:

            if x[1] in kegg_gene_to_path_codes:
                path2 = kegg_gene_to_path_codes[x[1]]

            if x[1] in reactome_gene_to_path_codes:
                path2.extend(reactome_gene_to_path_codes[x[1]])

        total = list(set(path1).union(path2))
        common = list(set(path1).intersection(path2))

        vqm = np.sqrt((len(path1) ** 2 + len(path2) ** 2) / 2)
        data[0][0] = vqm

        if len(total) == 0:
            data[0][1] = 0.
        else:
            data[0][1] = float(len(common)) / len(total)

        # HPO
        hpo1 = []
        hpo2 = []
        if x[0] in hpo_gene_to_code:
            hpo1 = hpo_gene_to_code[x[0]]
        if x[1] in hpo_gene_to_code:
            hpo2 = hpo_gene_to_code[x[1]]
        total = list(set(hpo1).union(hpo2))
        common = list(set(hpo1).intersection(hpo2))
        vqm = np.sqrt((len(hpo1) ** 2 + len(hpo2) ** 2) / 2)

        data[0][2] = vqm

        if len(total) == 0:
            data[0][3] = 0.
        else:
            data[0][3] = float(len(common)) / len(total)

        # PPI Network
        dist = []
        neighbors1 = []
        neighbors2 = []
        if x[0] in dists_ppi:
            neighbors1 = [p for p in nx.all_neighbors(G_ppi, x[0]) if p != x[0]]
            if x[1] in dists_ppi[x[0]]:
                dist.append(dists_ppi[x[0]][x[1]])
        if x[1] in dists_ppi:
            neighbors2 = [p for p in nx.all_neighbors(G_ppi, x[1]) if p != x[1]]
            if x[0] in dists_ppi[x[1]]:
                dist.append(dists_ppi[x[1]][x[0]])
        if dist != [] and min(dist) > 0:
            ppi_dist = 1 / float(min(dist))
        else:
            ppi_dist = 0.

        total = list(set(neighbors1).union(neighbors2))
        common = list(set(neighbors1).intersection(neighbors2))
        vqm = np.sqrt((len(neighbors1) ** 2 + len(neighbors2) ** 2) / 2)

        data[0][4] = vqm

        if len(total) == 0:
            data[0][5] = 0.
        else:
            data[0][5] = float(len(common)) / len(total)

        # data[i][8] = len(common)
        data[0][6] = ppi_dist

        # PWY Network
        dist = []
        neighbors1 = []
        neighbors2 = []
        if x[0] in dists_pwy:
            neighbors1 = [p for p in nx.all_neighbors(G_pwy, x[0]) if p is not x[0]]
            if x[1] in dists_pwy[x[0]]:
                dist.append(dists_pwy[x[0]][x[1]])
        if x[1] in dists_pwy:
            neighbors2 = [p for p in nx.all_neighbors(G_pwy, x[1]) if p is not x[1]]
            if x[0] in dists_pwy[x[1]]:
                dist.append(dists_pwy[x[1]][x[0]])

        if dist != [] and min(dist) > 0:
            pwy_dist = 1 / float(min(dist))
        else:
            pwy_dist = 0.

        total = list(set(neighbors1).union(neighbors2))
        common = list(set(neighbors1).intersection(neighbors2))
        vqm = np.sqrt((len(neighbors1) ** 2 + len(neighbors2) ** 2) / 2)

        data[0][7] = vqm

        if len(total) == 0:
            data[0][8] = 0.
        else:
            data[0][8] = float(len(common)) / len(total)

        # data[i][12] = len(common)
        data[0][9] = pwy_dist

        # TXT Network
        dist = []
        neighbors1 = []
        neighbors2 = []
        if x[0] in dists_txt:
            neighbors1 = [p for p in nx.all_neighbors(G_txt, x[0]) if p is not x[0]]
            if x[1] in dists_txt[x[0]]:
                dist.append(dists_txt[x[0]][x[1]])
        if x[1] in dists_txt:
            neighbors2 = [p for p in nx.all_neighbors(G_txt, x[1]) if p is not x[1]]
            if x[0] in dists_txt[x[1]]:
                dist.append(dists_txt[x[1]][x[0]])

        if dist != [] and min(dist) > 0:
            txt_dist = 1 / float(min(dist))
        else:
            txt_dist = 0.

        total = list(set(neighbors1).union(neighbors2))
        common = list(set(neighbors1).intersection(neighbors2))
        vqm = np.sqrt((len(neighbors1) ** 2 + len(neighbors2) ** 2) / 2)

        data[0][10] = vqm

        if len(total) == 0:
            data[0][11] = 0.
        else:
            data[0][11] = float(len(common)) / len(total)

        # data[i][16] = len(common)
        data[0][12] = txt_dist

        # Co-expression

        rankcoex1 = []
        rankcoex2 = []
        coexvalue = 0.
        if x[0] in coexpress_dict:
            rankcoex1 = [c for c in coexpress_dict[x[0]] if coexpress_dict[x[0]][c] < 100]
            if x[1] in coexpress_dict[x[0]]:
                coexvalue = 1 / coexpress_dict[x[0]][x[1]]
        if x[1] in coexpress_dict:
            rankcoex2 = [c for c in coexpress_dict[x[1]] if coexpress_dict[x[1]][c] < 100]
            if x[0] in coexpress_dict[x[1]]:
                coexvalue = 1 / coexpress_dict[x[1]][x[0]]

        total = list(set(rankcoex1).union(rankcoex2))
        common = list(set(rankcoex1).intersection(rankcoex2))
        vqm = np.sqrt((len(rankcoex1) ** 2 + len(rankcoex2) ** 2) / 2)

        data[0][13] = vqm

        if len(total) == 0:
            data[0][14] = 0.
        else:
            data[0][14] = float(len(common)) / len(total)

        # data[i][20] = len(common)
        data[0][15] = coexvalue

        # Lof

        if x[0] in lof_dict:
            v1 = lof_dict[x[0]]
        else:
            v1 = 0

        if x[1] in lof_dict:
            v2 = lof_dict[x[1]]
        else:
            v2 = 0

        vqm = np.sqrt((v1 ** 2 + v2 ** 2) / 2)
        data[0][16] = vqm

        # Happloinsufficiency Analysis

        if x[0] in hap_insuf_dict:
            v1 = hap_insuf_dict[x[0]]
        else:
            v1 = 0

        if x[1] in hap_insuf_dict:
            v2 = hap_insuf_dict[x[1]]
        else:
            v2 = 0

        vqm = np.sqrt((v1 ** 2 + v2 ** 2) / 2)
        data[0][17] = vqm

        # Protein Age

        if x[0] in protein_age_dict:
            v1 = protein_age_dict[x[0]]
        else:
            v1 = 0

        if x[1] in protein_age_dict:
            v2 = protein_age_dict[x[1]]
        else:
            v2 = 0

        vqm = np.sqrt((v1 ** 2 + v2 ** 2) / 2)
        data[0][18] = vqm

        # dN/DS

        if x[0] in dNdS_avg_dict:
            v1 = dNdS_avg_dict[x[0]]
        else:
            v1 = 0

        if x[1] in dNdS_avg_dict:
            v2 = dNdS_avg_dict[x[1]]
        else:
            v2 = 0

        vqm = np.sqrt((v1 ** 2 + v2 ** 2) / 2)
        data[0][19] = vqm

        # Gene Essentiality

        if x[0] in gene_ess_dict:
            v1 = np.mean(gene_ess_dict[x[0]])
        else:
            v1 = 0.
        if x[1] in gene_ess_dict:
            v2 = np.mean(gene_ess_dict[x[1]])
        else:
            v2 = 0.

        vqm = np.sqrt((v1 ** 2 + v2 ** 2) / 2)
        data[0][20] = vqm

    df = pd.DataFrame(all_data, index=new_list_pairs, columns=[
        # Pathways
        '#ofpathways',  # 0
        'common_pathways',  # 1
        # Phenotypes
        '#ofphenotypes',  # 2
        'common_phenotypes',  # 3
        # PPI network
        '#ofNeighborsPPI',  # 4
        '#Common_PPI_Neighbors',  # 5
        'PPI_network_dist',  # 6
        # PWY network
        '#ofNeighborsPWY',  # 7
        '#common_PWY_neighbors',  # 8
        'PWY_network_dist',  # 9
        # Txt network
        '#ofNeighborsTxt',  # 10
        '#Common_Txt_Neighbors',  # 11
        'Txt_network_dist',  # 12
        # Co-expression
        '#ofHighlyCoexpressed',  # 13
        '#Common_coexpressed',  # 14
        'Co-expression_coefficient',  # 15
        # LoF
        'LoFintolerance',  # 16
        # Haploinsuffiency
        'Haploinsufficiency',  # 17
        # Protein Age
        'protein_Age',  # 18
        # dN/dS
        'dN/dS',  # 19
        # Gene Essentiality
        'Essentiality',  # 20

    ])

    return df


if __name__ == '__main__':

    digepred_res_df = get_features(pairs)  # get feature values and save in a pandas DiGePred results df.
    p = clfs[model].predict_proba(digepred_res_df)[:, 1]  # get predictions based on DiGePred model specified.
    digepred_res_df[model] = p  # add column to DiGePred result df.

    digepred_res_df.to_csv('~/DiGePred/results/{path}_{project_name}.csv'.format(project_name=project_name),
                    sep=',', header=True, index=False)  # save feature values and predictions as DiGePred results CSV.

