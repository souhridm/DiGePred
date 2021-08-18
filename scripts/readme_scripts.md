Scripts to train models and test DiGePred performance.

DiGePred_train_all_sets.py: script to train the classifier and test performance during training. 

DiGePred_held_out_test_performance.py: script to get predictions on held-out test and measure performance using ROC and PR curves.

Get_DiGePred_scores_user_input_genes_or_pairs.py: script where user may provide an input file (.txt) containing either list of genes (-g) or list of gene pairs (-p), along with model of choice (-m) and job name (-n). It generates DiGePred results CSV with feature values and DiGePred predictions, based on model of choice. 

gene pairs should be specified as:
geneA,geneB (preferably in alphabetical order).

"unaffected-no-gene-overlap" is the default model, but user can select models among permuted, random, unaffected, random-no-gene-overlap and all-digenic-vs-unaffected.
