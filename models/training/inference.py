import os

import numpy as np
import torch
from torch.utils.data import DataLoader


def run_inference(model, dataset, device, out_dir, batch_size=128):
    """Run the trained model on every sample and save per-sample outputs.

    Outputs written to out_dir:
        latents.npy         — (n_samples, latent_dim)  mu vectors
        reconstructions.npy — (n_samples, n_bins)      log2(count+1) space
        sample_ids.npy      — (n_samples,)             sample ID strings
    """
    os.makedirs(out_dir, exist_ok=True)
    model.eval()

    # num_workers=0: required for MPS, safe everywhere for inference
    dl = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    all_mu    = []
    all_recon = []

    with torch.no_grad():
        for batch in dl:
            x     = batch.to(device)
            mu, _ = model.enc(x)
            recon = model.dec(mu)       # deterministic: use mu rather than sampled z
<<<<<<< HEAD
            all_mu.append(mu.cpu().numpy())
            all_recon.append(recon.cpu().numpy())
=======
            # Denormalize reconstructions: inverse of log2(count+1) is 2^x - 1
            recon_denorm = torch.pow(2, recon) - 1
            all_mu.append(mu.cpu().numpy())
            all_recon.append(recon_denorm.cpu().numpy())
>>>>>>> 0721fad3f32fa3b072c3211eef79f67cef0ed82a

    latents = np.concatenate(all_mu,    axis=0)   # (n_samples, latent_dim)
    recons  = np.concatenate(all_recon, axis=0)   # (n_samples, n_bins)

    np.save(os.path.join(out_dir, "latents.npy"),         latents)
    np.save(os.path.join(out_dir, "reconstructions.npy"), recons)
    np.save(os.path.join(out_dir, "sample_ids.npy"),      np.array(dataset.sample_ids))

    print(f"Saved latents      {latents.shape} → {out_dir}/latents.npy",          flush=True)
<<<<<<< HEAD
    print(f"Saved recons       {recons.shape}  → {out_dir}/reconstructions.npy",  flush=True)
=======
    print(f"Saved recons (denorm) {recons.shape}  → {out_dir}/reconstructions.npy",  flush=True)
>>>>>>> 0721fad3f32fa3b072c3211eef79f67cef0ed82a
    print(f"Saved sample_ids   ({len(dataset.sample_ids)},) → {out_dir}/sample_ids.npy", flush=True)
