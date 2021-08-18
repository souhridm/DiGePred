# DiGePred
Random Forest classifier designed to predict pairs of human genes capable to causing a digenic disease when carrying rare variants simultaneously. 

DiGePred has been trained using digenic pairs from DIDA (Gazzo et al., 2016; http://dida.ibsquare.be) and non-digenic pairs from various sources.

  unaffected: non-digenic pairs derived from unaffected relatives of individuals with rare undiagnosed diseases, as part of the Undiagnosed Diseases Network.
  permuted: non-digenic pairs derived from permutations of genes from digenic pairs in DIDA.
  random: non-digenic pairs generated randomly.
  matched: non-digenic pairs chosen to match the distribution of features of digenic pairs.
  unaffected no gene overlap: unaffected non-digenic pairs chosen such that no genes were common between the training and held-out test sets.  
  random no gene overlap: random non-digenic pairs chosen such that no genes were common between the training and held-out test sets.
 
The positive (digenic) and negative (non-digenic) gene pairs used to train and test DiGePred are provide in the flders "positives" and "negatives".

All trained DiGePred models are provided in the folder "models". "Unaffected-no-gene-overlap" was the final and best performing model. 

Scripts to train models and test DiGePred performance have been provided in "scripts".

  DiGePred_train_all_sets.py: script to train the classifier and test performance during training.
  DiGePred_held_out_test_performance.py: script to get predictions on held-out test and measure performance using ROC and PR curves.

DiGePred has been run on all human gene pairs, based on all genes from HGNC. The scores are available here:
  https://vanderbilt.box.com/s/n1nzdyj8i5fa55vultyq4xn6rsp792a7
  
  https://vanderbilt.box.com/s/459ethsqv339nqiarhm0j227jdjb0whq
  
  https://vanderbilt.box.com/s/acdqvjuihj3932c6msi5py82rvr5kam3
  
  https://vanderbilt.box.com/s/kb3vzubfxjcjtxt8x0y1vytu59x8r8no
  
A website is available where the user can access DiGePred scores for all human gene pairs (http://www.meilerlab.org/index.php/servers/show?s_id=28). 
