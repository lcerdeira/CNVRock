# Experiment 01 (manual)

Historical notes:
- Originally used 100bp bin size (so 230,000+ input size), so needed 50 GB VRAM. I then moved to 1000bp bins, before also filtering down to core genome. Models overfit so much so that real CNVs disappeared because the models "learned" the read count profiles containing CNVs, because I set patience to be way too high. 

Key setup highlights:
- patience 10, latent space 10
- drafting HMM classifier

Observations:
- latent variable 5 appears to correlate with the quality of the copy ratio signal. A high z5 (4-7) often means a really clean copy ratio hovering around 1. A low z5 (< -2) often show samples with 0 coverage or horrible sWGA amplification artefacts. 
- SPT38297 has CNV on chromosome 7. Potentially some kinases, transferases? And GCH1 triplication
- Poor recall on GCH1 amplifications

Things to try next:
- ought to do proper test-train-validation splits at some point
- Maybe need to start data augmentation and doing smart things with down-sampling on CN>1 samples