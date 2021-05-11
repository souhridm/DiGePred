import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
import datetime
from sklearn.metrics import *

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

clf_set_col = {'permuted': '#fd8d3c',
               'random': '#238b45',
               'matched': '#969696',
               'unaffected': '#6baed6',
               'all digenic vs unaffected': '#6baed6',
               'unaffected no gene overlap': '#08519c',
               'random no gene overlap': '#c7e9c0'
               }


models = {
    'unaffected': {'pos': pd.read_csv('~/digenic_pairs_held_out.csv')[sel_feats],
                   'neg': pd.read_csv('~/unaffected_non_digenic_pairs_held_out.csv')[sel_feats]
                   },
    'permuted': {'pos': pd.read_csv('~/digenic_pairs_held_out.csv')[sel_feats],
                 'neg': pd.read_csv('~/permuted_non_digenic_pairs_held_out.csv')[sel_feats]
                 },
    'random': {'pos': pd.read_csv('~/digenic_pairs_held_out.csv')[sel_feats],
               'neg': pd.read_csv('~/random_non_digenic_pairs_held_out.csv')[sel_feats]
               },
    'matched': {'pos': pd.read_csv('~/digenic_pairs_held_out.csv')[sel_feats],
                'neg': pd.read_csv('~/matched_non_digenic_pairs_held_out.csv')[sel_feats]
                },

    'unaffected no gene overlap': {'pos': pd.read_csv('~/digenic_pairs_no_overlap_held_out.csv')[sel_feats],
                                   'neg': pd.read_csv('~/unaffected no gene overlap_non_digenic_pairs_held_out.csv')[sel_feats]
                                   },
    'random no gene overlap': {'pos': pd.read_csv('~/digenic_pairs_no_overlap_held_out.csv')[sel_feats],
                               'neg': pd.read_csv('~/random no gene overlap_non_digenic_pairs_held_out.csv')[sel_feats]
                               },
            }

clfs = dict()
clfs['unaffected'] = pd.read_pickle(
    '~/DiGePred_unaffected_replace_ind_gene_feats_w QM_mean_Mar24_21.sav')
clfs['permuted'] = pd.read_pickle(
    '~/DiGePred_permuted_replace_ind_gene_feats_w QM_mean_Mar24_21.sav')
clfs['random'] = pd.read_pickle(
    '~/DiGePred_random_replace_ind_gene_feats_w QM_mean_Mar24_21.sav')
clfs['matched'] = pd.read_pickle(
    '~/DiGePred_matched_replace_ind_gene_feats_w QM_mean_Mar24_21.sav')

clfs['unaffected no gene overlap'] = pd.read_pickle(
    '~/Classifiers/DiGePred_unaffected no gene overlap_replace_ind_gene_feats_w QM_mean_Mar24_21.sav')
clfs['random no gene overlap'] = pd.read_pickle(
    '~/Classifiers/DiGePred_random no gene overlap_replace_ind_gene_feats_w QM_mean_Mar24_21.sav')


def get_held_out_validation_df(model):
    pos_validation_df = models[model]['pos']
    clf = clfs[model]
    neg_pairs = models[model]['neg']
    pos_neg_ratio = 75
    np.random.seed(seed=1607)

    min_num_neg = np.minimum(int(round(pos_neg_ratio * len(pos_validation_df))), len(neg_pairs))

    n_sets = np.minimum(len(neg_pairs.index) // min_num_neg, 10)

    indices = np.random.choice(len(neg_pairs.index), (n_sets, min_num_neg), replace=False)

    data_cols = {}
    for i, i_negs in enumerate(indices):

        data_cols[i] = {}
        non_digenic_df = neg_pairs.iloc[i_negs]

        X_validation = pd.concat([pos_validation_df, non_digenic_df], ignore_index=True)
        y_validation = np.asarray([1] * len(pos_validation_df) + [0] * len(non_digenic_df))

        preds = clf.predict_proba(X_validation)

        data_cols[i]['Test data'] = y_validation
        data_cols[i]['Pred probs'] = preds[:, 1]

    df = pd.DataFrame(index=range(len(indices)), columns=['Test data', 'Pred probs'])

    for i in data_cols:
        for p in data_cols[i]:
            df[p][i] = data_cols[i][p]

    return df


def get_pv(model, beta):
    pos_validation_df = models[model]['pos']
    np.random.seed(seed=1607)

    clf = clfs[model]
    neg_pairs = models[model]['neg']
    pos_neg_ratio = 75

    min_num_neg = np.minimum(int(round(pos_neg_ratio * len(pos_validation_df))), len(neg_pairs))

    n_sets = np.minimum(len(neg_pairs) // min_num_neg, 10)

    indices = np.random.choice(len(neg_pairs.index), (n_sets, min_num_neg), replace=False)

    data_cols = {}
    fbeta_scores = {}
    for i, i_negs in enumerate(indices):

        data_cols[i] = {}
        fbeta_scores[i] = {}
        non_digenic_df = neg_pairs.iloc[i_negs]

        X_validation = pd.concat([pos_validation_df, non_digenic_df], ignore_index=True)
        y_validation = np.asarray([1] * len(pos_validation_df) + [0] * len(non_digenic_df))

        preds = clf.predict_proba(X_validation)[:, 1]
        pred_w_labels = [(preds[i], y_validation[i]) for i in range(len(preds))]
        totp = len(pos_validation_df)
        totn = len(non_digenic_df)
        ntp = 0
        nfp = 0
        for (pred, label) in sorted(pred_w_labels, key=lambda tup: tup[0], reverse=True):
            if label == 1:
                ntp += 1
            else:
                nfp += 1

            nfn = totp - ntp
            ntn = totn - nfp

            tpr, fpr, precision, recall = get_tpr_fpr_precision_recall(ntp=ntp, nfp=nfp, nfn=nfn, ntn=ntn)
            if precision is not np.nan:
                fbeta_score = get_fx_score(precision=precision, recall=recall, beta=beta)
                if fbeta_score is not np.nan:
                    if round(pred, 3) not in fbeta_scores[i]:
                        fbeta_scores[i][round(pred, 3)] = []
                    fbeta_scores[i][round(pred, 3)].append(fbeta_score)

    score_rank = {}
    for i in fbeta_scores:
        for n, j in enumerate(sorted(fbeta_scores[i], key=lambda x: np.mean(fbeta_scores[i][x]), reverse=True)):
            if j not in score_rank:
                score_rank[j] = []
            score_rank[j].append(n + 1)

    pv = sorted(score_rank, key=lambda x: np.mean(score_rank[x]))[0]

    return pv, fbeta_scores, score_rank


def get_tpr_fpr_precision_recall(ntp, nfp, ntn, nfn):
    tpr = float(ntp) / (ntp + nfn)
    fpr = float(nfp / (nfp + ntn))
    precision = [float(ntp) / (ntp + nfp) if (ntp + nfp) > 0 else np.nan][0]
    recall = tpr

    return tpr, fpr, precision, recall


def get_fx_score(precision, recall, beta):
    if (precision + recall) > 0:
        f = (1 + np.square(beta)) * precision * recall / ((np.square(beta) * precision) + recall)
    else:
        f = np.nan
    return f


def get_roc_pr_curve_data(prediction_df):
    preds = []
    labels = []
    for i_pair in prediction_df.index:
        pr = prediction_df.loc[i_pair]['Pred probs']
        ls = prediction_df.loc[i_pair]['Test data']

        preds.append(pr)
        labels.append(ls)

    # Calculating metrics for ROC PR curves

    overall_tprs = []
    overall_fprs = []
    overall_precisions = []
    overall_recalls = []
    overall_preds = []
    pred_value_dict = {}

    for i in range(len(preds)):

        pred_w_labels = [(preds[i][j], labels[i][j]) for j in range(len(preds[i]))]

        n_identified_pos = 0
        n_identified_neg = 0

        n_tot_pos = list(labels[i]).count(1)
        n_tot_neg = list(labels[i]).count(0)

        tp = 0.
        fp = 0.

        all_tprs = []
        all_fprs = []
        all_precisions = []
        all_recalls = []
        all_preds = []

        for (p, l) in sorted(pred_w_labels, key=lambda tup: tup[0], reverse=True):

            if p not in pred_value_dict:
                pred_value_dict[p] = {'fpr': [], 'tpr': [], 'precision': [], 'recall': []}

            if l == 1:
                tp += 1
            else:
                fp += 1

            tn = n_tot_neg - fp
            fn = n_tot_pos - tp

            tpr = tp / (tp + fn)
            fpr = fp / (fp + tn)

            precision = tp / (tp + fp)
            recall = tpr

            for i in range(len(all_precisions)):
                vi = all_precisions[i]

                if vi < precision:
                    all_precisions[i] = precision

            all_precisions.append(precision)
            all_tprs.append(tpr)
            all_fprs.append(fpr)

            all_recalls.append(recall)
            all_preds.append(p)

        for i, p in enumerate(all_preds):

            pred_value_dict[p]['fpr'].append(all_fprs[i])
            pred_value_dict[p]['tpr'].append(all_tprs[i])
            pred_value_dict[p]['precision'].append(all_precisions[i])
            pred_value_dict[p]['recall'].append(all_recalls[i])

        overall_fprs.append(all_fprs)
        overall_tprs.append(all_tprs)
        overall_precisions.append(all_precisions)
        overall_recalls.append(all_recalls)
        overall_preds.append(all_preds)

    return overall_tprs, overall_fprs, overall_precisions, overall_recalls, overall_preds, pred_value_dict


def plot_roc_pr_curves(neg_sets):
    mpl.rc('font', family='Arial')
    fig, ax = plt.subplots(1, 2, figsize=(10, 7))

    ax[0].set_aspect('equal', 'box')
    ax[1].set_aspect('equal', 'box')

    leg_roc_handles = []
    leg_roc_labels = []

    leg_pr_handles = []
    leg_pr_labels = []

    for m in neg_sets:
        predicted_df = get_held_out_validation_df(model=m)

        overall_tprs, overall_fprs, overall_precisions, overall_recalls, overall_preds, pred_value_dict = \
            get_roc_pr_curve_data(prediction_df=predicted_df)

        mean_fpr = np.mean(overall_fprs, axis=0)
        mean_tpr = np.mean(overall_tprs, axis=0)
        mean_tpr[0] = 0

        std_tpr = np.std(overall_tprs, axis=0)

        roc_auc = auc(mean_fpr, mean_tpr)

        roc_line, = ax[0].step(mean_fpr, mean_tpr,
                               color=clf_set_col[m],
                               alpha=0.6,
                               linewidth=2,
                               linestyle='-')

        std_roc_auc = np.mean([auc(mean_fpr, np.minimum(mean_tpr + std_tpr, 1)) - roc_auc,
                               roc_auc - auc(mean_fpr, np.maximum(mean_tpr - std_tpr, 0))])

        leg_roc_handles.append((roc_line))
        leg_roc_labels.append('{ns} ({m})'
                              .format(m=round(roc_auc, 3), s=round(std_roc_auc, 3), ns=m))

        mean_recall = np.mean(overall_recalls, axis=0)
        mean_precision = np.mean(overall_precisions, axis=0)
        mean_recall[0] = 0

        std_precision = np.std(overall_precisions, axis=0)

        pr_auc = auc(mean_recall, mean_precision)

        pr_line, = ax[1].step(mean_recall, mean_precision,
                              color=clf_set_col[m],
                              alpha=0.6,
                              linewidth=2,
                              linestyle='-')

        std_pr_auc = np.mean([auc(mean_recall, np.minimum(mean_precision + std_precision, 1)) - pr_auc,
                              pr_auc - auc(mean_recall, np.maximum(mean_precision - std_precision, 0))])

        leg_pr_handles.append((pr_line))
        leg_pr_labels.append('{ns} ({m})'
                             .format(m=round(pr_auc, 3), s=round(std_pr_auc, 3), ns=m))

    c, = ax[0].plot([0, 1], [0, 1], linestyle='--', lw=1, color='k', alpha=0.6)
    leg_roc_handles.append(c)
    leg_roc_labels.append('chance ({})'.format(0.5))

    leg1 = ax[0].legend(leg_roc_handles, leg_roc_labels, loc='lower right', fontsize=8,
                        bbox_to_anchor=(1, 0.05), frameon=True, framealpha=0.1,
                        )
    leg1.set_title('held out gene pairs (mean AUC)', prop={'size': 10})
    ax[0].add_artist(leg1)

    c, = ax[1].plot([0, 1], [mean_precision[-1], mean_precision[-1]],
                    linestyle='--', lw=1, color='k', alpha=0.6)
    leg_pr_handles.append(c)
    leg_pr_labels.append('chance ({})'.format(0.0132))

    # legfp = ax[1].legend(legf_pr_handles, legf_pr_labels, loc='lower left', fontsize=12,
    #                      bbox_to_anchor=(0.5, 0.5, 0.5, 0.5))
    # legfp.set_title('Thresholds', prop={'size': 14})
    # ax[1].add_artist(legfp)

    leg2 = ax[1].legend(leg_pr_handles, leg_pr_labels, loc='lower left', fontsize=8,
                        bbox_to_anchor=(0, 0.05), frameon=True, framealpha=0.1,
                        )

    leg2.set_title('held out gene pairs (mean AUC)', prop={'size': 10})
    ax[1].add_artist(leg2)

    ax[0].set_xlim([-0.02, 1.02])
    ax[0].set_ylim([-0.02, 1.02])
    ax[0].set_xticks(np.linspace(0, 1, 6))
    ax[0].set_yticks(np.linspace(0, 1, 6))
    ax[0].set_xticklabels([round(v, 1) for v in np.linspace(0, 1, 6)], fontsize=10)
    ax[0].set_yticklabels([round(v, 1) for v in np.linspace(0, 1, 6)], fontsize=10)
    ax[0].set_xlabel('False Positive Rate', fontsize=12)
    ax[0].set_ylabel('True Positive Rate', fontsize=12)
    ax[0].set_title('ROC Curve', fontsize=12)
    ax[0].grid(b=True, which='major', color='#bdbdbd', linestyle='--', lw=0.7, alpha=0.5)

    ax[1].set_xlim([-0.02, 1.02])
    ax[1].set_ylim([-0.02, 1.02])
    ax[1].set_xticks(np.linspace(0, 1, 6))
    ax[1].set_yticks(np.linspace(0, 1, 6))
    ax[1].set_xticklabels([round(v, 1) for v in np.linspace(0, 1, 6)], fontsize=10)
    ax[1].set_yticklabels([round(v, 1) for v in np.linspace(0, 1, 6)], fontsize=10)
    ax[1].set_xlabel('Recall', fontsize=12)
    ax[1].set_ylabel('Precision', fontsize=12)
    ax[1].set_title('PR Curve', fontsize=12)
    ax[1].grid(b=True, which='major', color='#bdbdbd', linestyle='--', lw=0.7, alpha=0.5)

    fig.savefig('~/DiGePred_test_performance_{month}{day}_{year}.pdf'
                .format(month=month, day=day, year=year), bbox_inches='tight')
