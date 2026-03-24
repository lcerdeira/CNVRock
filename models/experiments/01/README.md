# Experiment 01

Model massively overfit. sWGA copy ratios are near 1, suggesting the decoder has fit over lots of sWGA profiles, so much so that real CNVs have disappeared. Likely need to restrict bottleneck. And let's only train on core regions, as these regions are extremely noisy anyway. Best to reduce model size. 

Weights are 41 MB. Reconstructions are 4.6 GB. 