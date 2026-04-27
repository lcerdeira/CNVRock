"""
06_conv_vae — genome-agnostic 1-D Convolutional VAE.
=====================================================

Identical architecture to 05_conv_vae (residual skip connections in encoder)
except that N_BINS_RAW is no longer a module-level constant — it is passed
as ``n_bins_raw`` to ConvVAE.__init__. This makes the architecture reusable
for any organism without changing code.

For P. falciparum use, pass n_bins_raw=20814 (or omit; defaults to 20814).
For K. pneumoniae HS11286, pass n_bins_raw to match the core genome bin count
extracted from your readcounts NPY (typically 4000–5000 bins).

The dataset automatically pads input to the next multiple of 32 (to satisfy
the 5 stride-2 layers), so the model's n_bins_padded must match. The dataset
derives this the same way: ``math.ceil(n_bins_raw / 32) * 32``.

Architecture (unchanged from 05):
  Encoder: 5× ResConvBlock, channels 1→32→64→128→256→256,
           main: Conv1d(k=7, stride=2) → BN → ReLU → Dropout(p=0.30),
           shortcut: Conv1d(k=1, stride=2) → BN.
  Decoder: 5× ConvTranspose1d blocks, channels 256→256→128→64→32→1,
           no dropout.
"""

import math
import torch
import torch.nn as nn


# ---------------------------------------------------------------------------
# Encoder block with residual skip connection (unchanged from 05_conv_vae)
# ---------------------------------------------------------------------------

class ResConvBlock(nn.Module):
    """Stride-2 conv block with a 1×1 strided shortcut to preserve spatial detail."""

    def __init__(self, in_ch: int, out_ch: int, stride: int = 2):
        super().__init__()
        self.main = nn.Sequential(
            nn.Conv1d(in_ch, out_ch, kernel_size=7, stride=stride, padding=3),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
            nn.Dropout(p=0.30),
        )
        self.shortcut = nn.Sequential(
            nn.Conv1d(in_ch, out_ch, kernel_size=1, stride=stride),
            nn.BatchNorm1d(out_ch),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.main(x) + self.shortcut(x)


# ---------------------------------------------------------------------------
# Encoder / Decoder
# ---------------------------------------------------------------------------

class ConvEncoder(nn.Module):
    def __init__(self, latent_dim: int, n_bins_padded: int):
        super().__init__()
        self.blocks = nn.Sequential(
            ResConvBlock(1,   32,  stride=2),
            ResConvBlock(32,  64,  stride=2),
            ResConvBlock(64,  128, stride=2),
            ResConvBlock(128, 256, stride=2),
            ResConvBlock(256, 256, stride=2),
        )
        self.flat_dim = 256 * (n_bins_padded // 32)
        self.mu     = nn.Linear(self.flat_dim, latent_dim)
        self.logvar = nn.Linear(self.flat_dim, latent_dim)

    def forward(self, x: torch.Tensor):
        h      = self.blocks(x.unsqueeze(1)).flatten(1)
        mu     = self.mu(h)
        logvar = torch.clamp(self.logvar(h), -10, 10)
        return mu, logvar


class ConvDecoder(nn.Module):
    def __init__(self, latent_dim: int, n_bins_raw: int, n_bins_padded: int):
        super().__init__()
        self._n_bins_raw    = n_bins_raw
        self._n_bins_padded = n_bins_padded
        self.flat_dim = 256 * (n_bins_padded // 32)
        self.proj     = nn.Linear(latent_dim, self.flat_dim)
        self.deconv   = nn.Sequential(
            self._block(256, 256, stride=2),
            self._block(256, 128, stride=2),
            self._block(128,  64, stride=2),
            self._block( 64,  32, stride=2),
            self._block( 32,   1, stride=2),
        )

    @staticmethod
    def _block(in_ch: int, out_ch: int, stride: int) -> nn.Sequential:
        return nn.Sequential(
            nn.ConvTranspose1d(in_ch, out_ch, kernel_size=7, stride=stride,
                               padding=3, output_padding=1),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        h     = self.proj(z).view(z.size(0), 256, self._n_bins_padded // 32)
        recon = self.deconv(h).squeeze(1)
        return recon[:, :self._n_bins_raw]


# ---------------------------------------------------------------------------
# VAE
# ---------------------------------------------------------------------------

_PF_N_BINS_RAW = 20814  # P. falciparum default — keeps existing Pf configs working

class ConvVAE(nn.Module):
    def __init__(self, latent_dim: int, n_bins_raw: int = _PF_N_BINS_RAW):
        super().__init__()
        n_bins_padded = math.ceil(n_bins_raw / 32) * 32
        self.enc = ConvEncoder(latent_dim, n_bins_padded)
        self.dec = ConvDecoder(latent_dim, n_bins_raw, n_bins_padded)

    def forward(self, x: torch.Tensor) -> dict:
        mu, logvar = self.enc(x)
        z          = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return {"recon": self.dec(z), "z": (mu, logvar)}
