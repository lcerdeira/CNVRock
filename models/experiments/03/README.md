# Experiment 03

Key setup highlights:
- patience 10, latent space 10

Observations:
- latent variable 5 appears to correlate with the quality of the copy ratio signal. A high z5 (4-7) often means a really clean copy ratio hovering around 1. A low z5 (< -2) often show samples with 0 coverage or horrible sWGA amplification artefacts. 

Things to try next:
- ought to do proper test-train-validation splits at some point
- it feels like I'm ready to build a classifier using an HMM and potentially using the outputs of the classifier to do something interesting with model training or data augmentation