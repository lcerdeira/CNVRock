import torch
import torch.nn.functional as F

from architectures import N_BINS_RAW


def compute_loss(x, outputs, beta):
    recon      = outputs["recon"]
    mu, logvar = outputs["z"]

    recon_loss = F.mse_loss(recon, x[:, :N_BINS_RAW], reduction="sum") / x.size(0)
    kl         = (-0.5 * (1 + logvar - mu.pow(2) - logvar.exp())).sum(1).mean()

    return recon_loss + beta * kl, {"recon": recon_loss, "kl": kl}


def train_vae(model, dataloader, optimiser,
              epochs, max_beta, warmup_epochs,
              patience, device, model_save_path=None):

    model.to(device)

    best_recon = float("inf")
    no_improve = 0

    try:
        for epoch in range(epochs):
            model.train()
            beta                     = max_beta * min(1.0, epoch / max(warmup_epochs, 1))
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
