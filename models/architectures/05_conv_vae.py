"""
05_conv_vae — 1-D Convolutional VAE with residual skip connections in encoder
==============================================================================

Architecture
------------
Identical to 04_conv_vae except that each stride-2 encoder block now has a
residual skip connection that adds a 1×1-strided shortcut to the main path:

    out = conv_bn_relu_dropout(x) + conv1x1_bn_stride(x)

Decoder is unchanged from 04_conv_vae.

Why skip connections
--------------------
Exp 12 (04_conv_vae, p=0.30, self_transition=0.80): PM2_PM3 FNR dropped from
0.93 to 0.53 — higher dropout strengthened anomaly signal. Exp 13 (same
checkpoint, self_transition=0.65): PM2_PM3 FNR fell further to 0.34, and FN
p50 = 1.26 >> 1.0 confirms real signal is still present.

However, GCH1 FNR is stuck at 0.11 across both exp 12 and 13 — completely
unresponsive to self_transition tuning (0.80 → 0.65). GCH1 FN p50 = 1.23
(unchanged between experiments) rules out HMM inertia as the GCH1 bottleneck.
HMM transitions at exp 13 p50=101 indicate the HMM is already quite jumpy;
further self_transition reduction risks over-segmentation.

The hypothesis: stride-2 convolutions discard fine spatial resolution that is
critical for detecting short (1–2 bin) amplification runs. PM2_PM3 and GCH1
have the shortest typical amplification spans; their anomaly signal gets
smeared across the stride-2 downsampling chain, making the per-bin CRR signal
inconsistent. Skip connections propagate the fine-scale input structure into
every level of the encoder, producing latent codes that are more sensitive to
spatially localised amplitude spikes. At inference, this yields a larger and
more spatially precise reconstruction error at the amplified bins, giving the
HMM a cleaner surface to commit to.

HMM self_transition is kept at 0.65 (the exp 13 value). Because emissions are
always re-fitted from the new checkpoint, keeping self_transition constant
isolates the effect of the architecture change.

Input / output: unchanged from 04_conv_vae (20814 bins, latent_dim configurable).
Encoder:  5× ResConvBlock, channels 1→32→64→128→256→256,
          main path: Conv1d(k=7, stride=2) → BN → ReLU → Dropout(p=0.30),
          shortcut:  Conv1d(k=1, stride=2) → BN.
Decoder:  5× ConvTranspose1d blocks — identical to 04_conv_vae, no dropout.
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
# Encoder block with residual skip connection
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
        # 1×1 conv to match output spatial size and channel count for the residual add
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
    def __init__(self, latent_dim: int):
        super().__init__()
        self.blocks = nn.Sequential(
            ResConvBlock(1,   32,  stride=2),
            ResConvBlock(32,  64,  stride=2),
            ResConvBlock(64,  128, stride=2),
            ResConvBlock(128, 256, stride=2),
            ResConvBlock(256, 256, stride=2),
        )
        self.flat_dim = 256 * (N_BINS_PADDED // 32)
        self.mu     = nn.Linear(self.flat_dim, latent_dim)
        self.logvar = nn.Linear(self.flat_dim, latent_dim)

    def forward(self, x: torch.Tensor):
        h      = self.blocks(x.unsqueeze(1)).flatten(1)
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
    def _block(in_ch: int, out_ch: int, stride: int) -> nn.Sequential:
        return nn.Sequential(
            nn.ConvTranspose1d(in_ch, out_ch, kernel_size=7, stride=stride,
                               padding=3, output_padding=1),
            nn.BatchNorm1d(out_ch),
            nn.ReLU(),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
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

    def forward(self, x: torch.Tensor) -> dict:
        mu, logvar = self.enc(x)
        z          = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return {"recon": self.dec(z), "z": (mu, logvar)}
