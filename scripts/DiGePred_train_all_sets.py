from sklearn.metrics import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedShuffleSplit
import numpy as np
import pandas as pd

import datetime
import random

now = datetime.datetime.now()
month = str(now.strftime("%b"))
day = str(now.strftime("%d"))
year = str(now.strftime("%y"))

sel_feats = ['common_pathways',
           'common_phenotypes',
           'Co-expression_coefficient',
           'PPI_network_dist',
           'PWY_network_dist',
           'Txt_network_dist',
           'LoFintolerance_combined',
           'Haploinsufficiency_combined',
           'Age_combined',
           'dN/dS_combined',
           'Essentiality_combined',
           '#ofpathways_combined',
           '#ofPhenotypeCodes_combined',
           '#ofNeighborsPPI_combined',
           '#ofNeighborsPWY_combined',
           '#ofNeighborsTxt_combined',
           '#ofHighlyCoexpressed_combined',
           '#Common_PPI_Neighbors',
           '#common_PWY_neighbors',
           '#Common_Txt_Neighbors',
           '#Common_coexpressed',
           ]

# Import dfs

digenic_training = pd.read_csv('~/')
digenic_training_no_overlap = pd.read_csv('~/')
unaffected_non_digenic_training = pd.read_csv('~/')
random_non_digenic_training = pd.read_csv('~/')
permuted_non_digenic_training = pd.read_csv('~/')
matched_non_digenic_training = pd.read_csv('~/')
unaffected_no_gene_overlap_non_digenic_training = pd.read_csv('~/')
random_no_gene_overlap_non_digenic_training = pd.read_csv('~/')

# define models with correpsonding positive and negative sets

models = {'unaffected': {'pos': digenic_training,
                         'neg': unaffected_non_digenic_training},
          'permuted': {'pos': digenic_training,
                       'neg': permuted_non_digenic_training},
          'random': {'pos': digenic_training,
                     'neg': random_non_digenic_training},
          'matched': {'pos': digenic_training,
                      'neg': matched_non_digenic_training},
          'unaffected no gene overlap': {'pos': digenic_training_no_overlap,
                                         'neg': unaffected_no_gene_overlap_non_digenic_training},
          'random no gene overlap': {'pos': digenic_training_no_overlap,
                                     'neg': random_no_gene_overlap_non_digenic_training},
          }

roc_aucs = {}
pr_aucs = {}
f1_scores = {}
tprs = {}
fprs = {}
precisions = {}
recalls = {}
test_data = {}
pred_probs = {}

pos_neg_ratio = 75

# Looping through all models


for m in models:

    digenic_pairs_df = models[m]['pos']
    digenic_pairs_list = list(digenic_pairs_df.index)

    non_digenic_pairs_df = models[m]['neg']
    non_digenic_pairs_list = list(non_digenic_pairs_df.index)

    roc_aucs[m] = []
    pr_aucs[m] = []
    f1_scores[m] = []

    tprs[m] = []
    fprs[m] = []
    precisions[m] = []
    recalls[m] = []

    test_data[m] = []
    pred_probs[m] = []

    # Making full set by adding digenic and non-digenic sets

    full_set_X = pd.concat([digenic_pairs_df, non_digenic_pairs_df], ignore_index=True)
    full_set_y = np.asarray([1] * len(digenic_pairs_list) + [0] * len(non_digenic_pairs_list))

    sss_splitter = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=xxxx)

    clf = RandomForestClassifier(n_jobs=1, n_estimators=500, max_depth=15)

    for train, test in sss_splitter.split(full_set_X, full_set_y):
        X_train, X_test, y_train, y_test = full_set_X.iloc[train], \
                                           full_set_X.iloc[test], \
                                           full_set_y[train], \
                                           full_set_y[test]

        # Train classifier using inner train

        clf.fit(X_train, y_train)
        preds = clf.predict_proba(X_test)
        predictions = clf.predict(X_test)

        # Evaluate performance using ROC, PR AUCs and F1 scores

        fpr, tpr, thresholds = roc_curve(y_test, preds[:, 1])
        roc_auc = auc(fpr, tpr)

        p, r, _ = precision_recall_curve(y_test, preds[:, 1])
        ap = average_precision_score(y_test, preds[:, 1])
        # f1 = (2 * p * r) / (p + r)

        f1 = f1_score(y_test, predictions)

        roc_aucs[m].append(roc_auc)
        pr_aucs[m].append(ap)
        f1_scores[m].append(f1)

        tprs[m].append(tpr)
        fprs[m].append(fpr)
        precisions[m].append(p)
        recalls[m].append(r)

        test_data[m].append(y_test)
        pred_probs[m].append(preds[:, 1])

data_cols = {
    'ROC_AUCs': roc_aucs,
    'PR_AUCs': pr_aucs,
    'F1_scores': f1_scores,
    'TPRs': tprs,
    'FPRs': fprs,
    'Precisions': precisions,
    'Recalls': recalls,
    'Test data': test_data,
    'Pred probs': pred_probs
}

df = pd.DataFrame(index=list(models), columns=data_cols.keys())

for m in models:
    for d in data_cols:
        df[d][m] = data_cols[d][m]

df.to_pickle('~/DiGePred_training_performance_{month}{day}_{year}.pkl'.format(month=month, day=day, year=year))

