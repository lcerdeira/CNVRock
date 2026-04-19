"""
04_conv_vae — 1-D Convolutional VAE with increased encoder dropout
==================================================================

Architecture
------------
Identical to 02_conv_vae except that encoder conv blocks use
Dropout(p=0.30) instead of Dropout(p=0.20).

Why higher dropout
------------------
02_conv_vae (p=0.20) is the current best baseline: CRT FNR=0.02,
MDR1 FNR=0.03, GCH1 FNR=0.07, but PM2_PM3 FNR=0.93–0.94 persists
across exp 09/10 despite HMM self_transition tuning (0.80→0.75).
FN CRR p50=1.43 >> 1.0 confirms signal IS present — the bottleneck
is either signal magnitude or consistency for the very short
(estimated 1–2 bin) PM2_PM3 amplification runs.

Exp 11 tested p=0.10 (lower dropout) to improve per-bin consistency;
it caused severe regressions across all genes. The confound: changing
dropout changes the reconstruction error distribution, invalidating
the HMM emission parameters calibrated on a different distribution.
p=0.10 + self_transition=0.80 simply wasn't a fair comparison.

p=0.30 tests the opposite direction. More aggressive regularisation
forces the encoder to learn a tighter, more generalised representation
of normal coverage. At inference (dropout disabled), the reconstruction
becomes a "cleaner" normal-mode estimate — CNV samples produce larger
absolute reconstruction errors and therefore stronger anomaly signal.
The hypothesis is that stronger anomaly signal boosts the signal-to-noise
ratio for PM2_PM3's short runs enough for the HMM to commit.

HMM self_transition is reset to 0.80 (exp 09 baseline) to isolate
the dropout effect — a different dropout will produce a different
emission distribution, so p=0.30 needs its own self_transition
evaluation before any further HMM tuning is meaningful.

Input / output: unchanged from 02_conv_vae.
Encoder:  5× Conv1d blocks, channels 1→32→64→128→256→256, each followed by
          BatchNorm1d → ReLU → Dropout(p=0.30).
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
            nn.Dropout(p=0.30),  # increased from 0.20 to sharpen anomaly signal for PM2_PM3
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
