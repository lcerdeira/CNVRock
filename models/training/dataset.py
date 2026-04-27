import math
import numpy as np
import torch
import torch.nn.functional as F
import os
from torch.utils.data import Dataset


class ReadCountDataset(Dataset):
    def __init__(self, npy_dir: str, normalise: bool = True):
        # npy_dir should contain counts.npy and sample_ids.npy
        self.counts     = np.load(os.path.join(npy_dir, "counts.npy"))
        self.sample_ids = np.load(os.path.join(npy_dir, "sample_ids.npy"), allow_pickle=True)
        self.normalise  = normalise
        n_bins_raw = self.counts.shape[1]
        # Must be a multiple of 32 so 5 stride-2 encoder layers divide evenly
        self._n_bins_padded = math.ceil(n_bins_raw / 32) * 32

    def __len__(self):
        return len(self.sample_ids)

    def __getitem__(self, idx):
        counts = self.counts[idx, :].astype(np.float32)
        if self.normalise:
            counts = np.log2(counts + 1)
        t = torch.from_numpy(counts)
        t = F.pad(t, (0, self._n_bins_padded - t.shape[0]))
        return t
