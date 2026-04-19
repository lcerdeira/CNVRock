"""
03_conv_vae — 1-D Convolutional VAE with reduced encoder dropout
================================================================

Architecture
------------
Identical to 02_conv_vae except that encoder conv blocks use
Dropout(p=0.10) instead of Dropout(p=0.20).

Why lower dropout
-----------------
02_conv_vae (p=0.20) achieved excellent GCH1 recovery but PM2_PM3 FNR
remained catastrophically high (0.93–0.94 across exp 09/10).  The FN CRR
p50=1.43 >> 1.0 confirms signal IS present in the reconstruction, but the
HMM cannot commit to the very short PM2_PM3 amplification runs (estimated
1–2 bins).  HMM self_transition tuning (0.80→0.75→0.70) had negligible
effect, suggesting the bottleneck is signal consistency rather than HMM
threshold.

Dropout at p=0.20 introduces stochastic noise that causes PM2_PM3
amplification runs to appear intermittent rather than sustained across
adjacent bins — the HMM therefore never sees enough consecutive elevated
bins to trigger a state transition.  Lowering to p=0.10 reduces this
noise, making the elevated CRR signal in PM2_PM3 amplifications appear
more consistently across bins so the HMM can commit.

We keep encoder dropout (not zero) to preserve the distributed
representation benefit that recovered GCH1 in exp 09.  p=0.10 is a
targeted reduction, not a removal.

Input / output: unchanged from 02_conv_vae.
Encoder:  5× Conv1d blocks, channels 1→32→64→128→256→256, each followed by
          BatchNorm1d → ReLU → Dropout(p=0.10).
Decoder:  5× ConvTranspose1d blocks — same as 02_conv_vae, no dropout.
"""

import math
import torch
import torch.nn as nn


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BINS_RAW    = 20814
N_BINS_PADDED = math.ceil(N_BINS_RAW / 32) * 32   # → 20832


# ---------------------------------------------------------------------------
# Encoder / Decoder blocks
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
            nn.Dropout(p=0.10),  # reduced from 0.20 to improve PM2_PM3 signal consistency
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


# ---------------------------------------------------------------------------
# VAE
# ---------------------------------------------------------------------------

class ConvVAE(nn.Module):
    def __init__(self, latent_dim: int):
        super().__init__()
        self.enc = ConvEncoder(latent_dim)
        self.dec = ConvDecoder(latent_dim)

    def forward(self, x):
        mu, logvar = self.enc(x)
        z          = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return {"recon": self.dec(z), "z": (mu, logvar)}
