import torch
import torch.nn as nn


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

N_BINS_RAW    = 23340
N_BINS_PADDED = 23360   # next multiple of 32


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
