import numpy as np
import torch
import torch.nn.functional as F
import zarr
from torch.utils.data import Dataset

from architectures.conv_vae import N_BINS_RAW, N_BINS_PADDED


class ReadCountDataset(Dataset):
    def __init__(self, store_path: str, normalise: bool = True):
        store           = zarr.open(store_path, mode="r")
        self.counts     = store["counts"]
        self.sample_ids = store.attrs["sample_ids"]
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
