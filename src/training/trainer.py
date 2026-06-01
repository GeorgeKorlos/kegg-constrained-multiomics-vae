import json
import torch
import hashlib
from pathlib import Path
from torch.optim import Adam
from torch.utils.data import TensorDataset, DataLoader


from src.model.vae import MultiOmicsVAE
from src.training.losses import (
    reconstruction_loss_trans,
    reconstruction_loss_meta,
    kl_loss,
    compute_eta,
)

F_TRANS = 17384
F_META = 225
LATENT_DIM = 128


def make_split(n, out_path, seed=0, ratios=(0.8, 0.1, 0.1)):
    out_path = Path(out_path)

    if out_path.exists():
        with open(out_path, "r") as f:
            return json.load(f)

    g = torch.Generator().manual_seed(seed)
    perm = torch.randperm(n, generator=g).tolist()

    n_val = int(n * ratios[1])
    n_test = int(n * ratios[2])
    n_train = n - n_val - n_test

    train = perm[:n_train]
    val = perm[n_train : n_train + n_val]
    test = perm[n_train + n_val :]

    split = {
        "seed": seed,
        "ratios": list(ratios),
        "n": n,
        "train": train,
        "val": val,
        "test": test,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)

    json_bytes = json.dumps(split, sort_keys=True, indent=2).encode("utf-8")

    with open(out_path, "wb") as f:
        f.write(json_bytes)

    sha256 = hashlib.sha256(json_bytes).hexdigest()

    with open(str(out_path) + ".sha256", "w") as f:
        f.write(sha256)

    return split


def guard2_block_variance(latents, partition, alpha, counters, p):
    dim_var = latents.var(dim=0, unbiased=False)
    block_var = torch.stack([dim_var[idxs].mean() for idxs in partition])

    global_mean = block_var.mean()

    threshold = alpha * global_mean
    flagged = block_var < threshold

    updated_counters = []
    for flag, count in zip(flagged.tolist(), counters):
        updated_counters.append(count + 1 if flag else 0)

    failed = any(count >= p for count in updated_counters)

    return failed, block_var, updated_counters


def guard3_r2(x_true, x_hat, eps=1e-8):
    ss_res = ((x_true - x_hat) ** 2).sum()
    ss_tot = ((x_true - x_true.mean(dim=0)) ** 2).sum()
    r2 = 1 - ss_res / (ss_tot + eps)
    return r2


def make_partition_per_dim(latent_dim=128):
    return [[i] for i in range(latent_dim)]


def evaluate(model, x_trans, x_meta, beta, partition, guard2_counters, alpha, p):
    model.eval()
    with torch.no_grad():
        out = model(x_trans, x_meta, sample=False)

        L_trans = reconstruction_loss_trans(out["x_trans_hat"], x_trans)
        L_meta = reconstruction_loss_meta(out["x_meta_hat"], x_meta)
        L_kl = kl_loss(out["mu_joint"], out["logvar_joint"])

        val_loss = L_trans + L_meta + beta * L_kl

        guard2_failed, block_var, guard2_counters = guard2_block_variance(
            latents=out["mu_joint"],
            partition=partition,
            alpha=alpha,
            counters=guard2_counters,
            p=p,
        )

        r2_trans = guard3_r2(x_trans, out["x_trans_hat"])
        r2_meta = guard3_r2(x_meta, out["x_meta_hat"])

        guard3_trans_failed = r2_trans <= 0
        guard3_meta_failed = r2_meta <= 0

    return {
        "val_loss": float(val_loss.item()),
        "r2_trans": float(r2_trans),
        "r2_meta": float(r2_meta),
        "guard2_failed": bool(guard2_failed),
        "guard3_trans_failed": bool(guard3_trans_failed),
        "guard3_meta_failed": bool(guard3_meta_failed),
        "block_var": block_var.cpu(),
        "guard2_counters": guard2_counters,
    }


def train(model, data, config):

    split_path = Path(config["split_path"])
    sha_path = Path(str(split_path) + ".sha256")

    if not split_path.exists():
        raise FileNotFoundError(f"Missing split artifact: {split_path}")

    if not sha_path.exists():
        raise FileNotFoundError(f"Missing split checksum: {sha_path}")

    with open(split_path, "r") as f:
        split = json.load(f)

    split_sha256 = sha_path.read_text().strip()

    beta = config.get("beta", 1.0)

    alpha = config["alpha"]
    p = config["p"]

    eta_clamp = config.get(
        "eta_clamp",
        (0.1, 10.0),
    )

    patience = config["patience"]
    min_delta = config["min_delta"]

    stop_on_guard3_fail = config.get(
        "stop_on_guard3_fail",
        True,
    )

    batch_size = config["batch_size"]
    lr = config["lr"]
    epochs = config["epochs"]

    trans = torch.as_tensor(
        data["trans"],
        dtype=torch.float32,
    )

    meta = torch.as_tensor(
        data["meta"],
        dtype=torch.float32,
    )

    assert trans.shape[1] == F_TRANS
    assert meta.shape[1] == F_META
    assert LATENT_DIM == 128

    train_idx = split["train"]
    val_idx = split["val"]

    train_trans = trans[train_idx]
    train_meta = meta[train_idx]

    val_trans = trans[val_idx]
    val_meta = meta[val_idx]

    train_dataset = TensorDataset(
        train_trans,
        train_meta,
    )

    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
    )

    optimizer = Adam(
        model.parameters(),
        lr=lr,
    )

    partition = make_partition_per_dim(LATENT_DIM)

    guard2_counters = [0] * len(partition)

    best_val = float("inf")
    patience_ctr = 0

    eta_log = []
    history = []

    Path("models").mkdir(
        parents=True,
        exist_ok=True,
    )

    for epoch in range(epochs):

        model.train()

        epoch_eta = []

        for x_trans, x_meta in train_loader:

            out = model(
                x_trans,
                x_meta,
                sample=True,
            )

            z = out["z"]

            L_trans = reconstruction_loss_trans(
                out["x_trans_hat"],
                x_trans,
            )

            L_meta = reconstruction_loss_meta(
                out["x_meta_hat"],
                x_meta,
            )

            L_kl = kl_loss(
                out["mu_joint"],
                out["logvar_joint"],
            )

            eta = compute_eta(
                L_trans,
                L_meta,
                z,
            )

            eta_value = float(eta.detach().item())

            eta_log.append(eta_value)
            epoch_eta.append(eta_value)

            eta_c = eta.clamp(
                min=eta_clamp[0],
                max=eta_clamp[1],
            )

            total = (1.0 / eta_c) * L_trans + L_meta + beta * L_kl

            optimizer.zero_grad()
            total.backward()
            optimizer.step()
        metrics = evaluate(
            model=model,
            x_trans=val_trans,
            x_meta=val_meta,
            beta=beta,
            partition=partition,
            guard2_counters=guard2_counters,
            alpha=alpha,
            p=p,
        )

        guard2_counters = metrics["guard2_counters"]

        eta_epoch_mean = sum(epoch_eta) / len(epoch_eta)

        history.append(
            {
                "epoch": epoch,
                **metrics,
                "eta_epoch_mean": eta_epoch_mean,
            }
        )

        if metrics["guard2_failed"]:
            print(f"[WARN] Guard 2 triggered " f"at epoch {epoch}")

        guard3_failed = metrics["guard3_trans_failed"] or metrics["guard3_meta_failed"]

        if guard3_failed and stop_on_guard3_fail:
            raise RuntimeError(
                "Guard 3 failure. "
                f"R²_trans={metrics['r2_trans']:.4f}, "
                f"R²_meta={metrics['r2_meta']:.4f}"
            )

        val_loss = metrics["val_loss"]

        if val_loss < (best_val - min_delta):

            checkpoint = {
                "state_dict": model.state_dict(),
                "split_sha256": split_sha256,
                "config": config,
                "epoch": epoch,
                "metrics": metrics,
            }

            torch.save(
                checkpoint,
                "models/baseline_vae.pt",
            )

            best_val = val_loss
            patience_ctr = 0

        else:
            patience_ctr += 1

            if patience_ctr >= patience:
                print(f"Early stopping at epoch {epoch}")
                break

    return {
        "history": history,
        "eta_log": eta_log,
        "best_val": best_val,
    }
