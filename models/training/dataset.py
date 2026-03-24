import numpy as np
import torch
import torch.nn.functional as F
import os
from torch.utils.data import Dataset

from architectures.conv_vae import N_BINS_RAW, N_BINS_PADDED


class ReadCountDataset(Dataset):
    def __init__(self, npy_dir: str, normalise: bool = True):
        # npy_dir should contain counts.npy and sample_ids.npy
        self.counts     = np.load(os.path.join(npy_dir, "counts.npy"))
        self.sample_ids = np.load(os.path.join(npy_dir, "sample_ids.npy"), allow_pickle=True)
        self.normalise  = normalise

    def __len__(self):
        return len(self.sample_ids)

    def __getitem__(self, idx):
        counts = self.counts[idx, :].astype(np.float32)
        if self.normalise:
            counts = np.log2(counts + 1)
        t = torch.from_numpy(counts)
        t = F.pad(t, (0, N_BINS_PADDED - N_BINS_RAW))
        return t
