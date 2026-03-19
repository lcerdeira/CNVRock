import os
import zarr
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data import Dataset, DataLoader


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BINS_RAW    = 233337
N_BINS_PADDED = 233344  # next multiple of 32


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class ConvEncoder(nn.Module):
    def __init__(self, latent_dim: int):
        super().__init__()
        self.conv = nn.Sequential(
            self._block(1,   32,  stride=2),
            self._block(32,  64,  stride=2),
            self._block(64,  128, stride=2),
            self._block(128, 256, stride=2),
            self._block(256, 256, stride=2),
        )
        self.flat_dim = 256 * (N_BINS_PADDED // 32)
        self.mu     = nn.Linear(self.flat_dim, latent_dim)
        self.logvar = nn.Linear(self.flat_dim, latent_dim)

    @staticmethod
    def _block(in_ch, out_ch, stride):
        return nn.Sequential(
            nn.Conv1d(in_ch, out_ch, kernel_size=7, stride=stride, padding=3),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
        )

    def forward(self, x):
        h      = self.conv(x.unsqueeze(1)).flatten(1)
        mu     = self.mu(h)
        logvar = torch.clamp(self.logvar(h), -10, 10)
        return mu, logvar


class ConvDecoder(nn.Module):
    def __init__(self, latent_dim: int):
        super().__init__()
        self.flat_dim = 256 * (N_BINS_PADDED // 32)
        self.proj     = nn.Linear(latent_dim, self.flat_dim)
        self.deconv   = nn.Sequential(
            self._block(256, 256, stride=2),
            self._block(256, 128, stride=2),
            self._block(128,  64, stride=2),
            self._block( 64,  32, stride=2),
            self._block( 32,   1, stride=2),
        )

    @staticmethod
    def _block(in_ch, out_ch, stride):
        return nn.Sequential(
            nn.ConvTranspose1d(in_ch, out_ch, kernel_size=7, stride=stride,
                               padding=3, output_padding=1),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
        )

    def forward(self, z):
        h     = self.proj(z).view(z.size(0), 256, N_BINS_PADDED // 32)
        recon = self.deconv(h).squeeze(1)
        return recon[:, :N_BINS_RAW]


class ConvVAE(nn.Module):
    def __init__(self, latent_dim: int):
        super().__init__()
        self.enc = ConvEncoder(latent_dim)
        self.dec = ConvDecoder(latent_dim)

    def forward(self, x):
        mu, logvar = self.enc(x)
        z          = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return {"recon": self.dec(z), "z": (mu, logvar)}


# ---------------------------------------------------------------------------
# Loss
# ---------------------------------------------------------------------------

def compute_loss(x, outputs, beta):
    recon      = outputs["recon"]
    mu, logvar = outputs["z"]

    recon_loss = F.mse_loss(recon, x[:, :N_BINS_RAW], reduction="sum") / x.size(0)
    kl         = (-0.5 * (1 + logvar - mu.pow(2) - logvar.exp())).sum(1).mean()

    return recon_loss + beta * kl, {"recon": recon_loss, "kl": kl}


# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train_vae(model, dataloader, optimiser,
              epochs, max_beta, warmup_epochs,
              patience, model_save_path="../data/conv_vae.pth"):

    assert torch.cuda.is_available(), "CUDA not available."
    device = torch.device("cuda")
    model.to(device)

    best_recon = float("inf")
    no_improve = 0

    try:
        for epoch in range(epochs):
            model.train()
            beta                     = max_beta * min(1.0, epoch / warmup_epochs)
            total, tot_recon, tot_kl = 0, 0, 0

            for i, batch in enumerate(dataloader):
                x   = batch.to(device)
                out = model(x)
                loss, det = compute_loss(x, out, beta)

                optimiser.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                optimiser.step()

                total     += loss.item()
                tot_recon += det["recon"].item()
                tot_kl    += det["kl"].item()

                if epoch == 0 or i % 50 == 0:
                    print(f"  epoch {epoch+1} | batch {i}/{len(dataloader)} | loss {loss.item():.4f}", flush=True)

            print(f"{epoch+1:04d} | loss {total:.2f} | recon {tot_recon:.2f} | kl {tot_kl:.2f} | beta {beta:.3f}", flush=True)

            if tot_recon < best_recon - 1e-3:
                best_recon = tot_recon
                no_improve = 0
                if model_save_path:
                    torch.save(model.state_dict(), model_save_path)
            else:
                no_improve += 1
                if no_improve >= patience:
                    print(f"Early stopping at epoch {epoch+1}.")
                    break

    except KeyboardInterrupt:
        print("Interrupted — returning model.")

    return model


# ---------------------------------------------------------------------------
# Setup + run
# ---------------------------------------------------------------------------

def setup_and_train(
    store_path      = "../data/Pf9-53973-samples-100bp.zarr",
    latent_dim      = 16,
    batch_size      = 128,
    epochs          = 5000,
    max_beta        = 1.0,
    warmup_epochs   = 400,
    patience        = 50,
    model_save_path = "../data/conv_vae.pth",
):
    lsf_jobid = os.getenv("LSB_JOBID")
    if lsf_jobid:
        print(f"LSF job: {lsf_jobid}\n{'='*50}")

    ds  = ReadCountDataset(store_path, normalise=True)
    dl  = DataLoader(ds, batch_size=batch_size, shuffle=True, num_workers=8)

    device    = torch.device("cuda")
    model     = ConvVAE(latent_dim=latent_dim).to(device)
    optimiser = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

    print(f"n_bins={N_BINS_RAW} (padded to {N_BINS_PADDED}) | latent_dim={latent_dim} | samples={len(ds)}", flush=True)
    print(model)

    return train_vae(
        model, dl, optimiser,
        epochs          = epochs,
        max_beta        = max_beta,
        warmup_epochs   = warmup_epochs,
        patience        = patience,
        model_save_path = model_save_path,
    )


# ---------------------------------------------------------------------------
# Inference helpers
# ---------------------------------------------------------------------------

def get_latent_embeddings(model, store_path, batch_size=128):
    """Returns (n_samples, latent_dim) mu vectors."""
    device = torch.device("cuda")
    model.eval()
    dl  = DataLoader(ReadCountDataset(store_path), batch_size=batch_size, shuffle=False)
    mus = []
    with torch.no_grad():
        for batch in dl:
            mu, _ = model.enc(batch.to(device))
            mus.append(mu.cpu().numpy())
    return np.concatenate(mus, axis=0)


def reconstruct(model, store_path, batch_size=128):
    """Returns reconstructions in raw count space."""
    device = torch.device("cuda")
    model.eval()
    dl     = DataLoader(ReadCountDataset(store_path), batch_size=batch_size, shuffle=False)
    recons = []
    with torch.no_grad():
        for batch in dl:
            recon = model(batch.to(device))["recon"]
            recons.append(recon.cpu().numpy())
    log2_recon = np.concatenate(recons, axis=0)
    return 2 ** log2_recon - 1


if __name__ == "__main__":
    setup_and_train()