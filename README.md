# DiGePred
Random Forest classifier designed to predict pairs of human genes capable to causing a digenic disease when carrying rare variants simultaneously. DiGePred has been trained using digenic pairs from DIDA and non-digenic pairs from unaffected relatives of individuals with rare undiagnosed disease. 
Training data and classifiers trained using varying non-digenic sets are available:
 digenic: digenic pairs donwloaded form DIDA (Gazzo et al., 2016; http://dida.ibsquare.be)
 unaffected: non-digenic pairs derived from unaffected relatives of individuals with rare undiagnosed diseases, as part of the Undiagnosed Diseases Network.
 permuted: non-digenic pairs derived from permutations of DIDA.
 random: non-digenic pairs generated randomly.
 matched: non-digenic pairs choen to match the distribution of features of digenic pairs.
 unaffected no gene overlap: unaffected non-digenic pairs chosen such that no genes wer ecommon between the training and held-out test sets.
 random no gene overlap: random non-digenic pairs chosen such that no genes wer ecommon between the training and held-out test sets.
 
Held-out test data also provided for all sets.
DiGePred_train_all_sets.py: script to train the classifier and test performance during training.
DiGePred_held_out_test_performance.py: script to get predictions on held-out test and measure performance using ROC and PR curves.
