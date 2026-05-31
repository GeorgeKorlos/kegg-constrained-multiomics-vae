import torch
import torch.nn as nn

LATENT_DIM = 128
LOGVAR_MIN, LOGVAR_MAX = -10.0, 10.0


class _MLPEncoder(nn.Module):

    def __init__(self, input_dim, hidden_dims, latent_dim=LATENT_DIM):
        super().__init__()
        layers = []
        prev = input_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.LayerNorm(h), nn.GELU()]
            prev = h
        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(prev, 2 * latent_dim)
        self.latent_dim = latent_dim

    def forward(self, x):
        mu, logvar = self.head(self.backbone(x)).chunk(2, dim=-1)
        logvar = torch.clamp(logvar, LOGVAR_MIN, LOGVAR_MAX)
        return (mu, logvar)


class TranscriptomicsEncoder(_MLPEncoder):
    """17384 -> 1024 -> 512 -> 2x128"""

    def __init__(self, input_dim=17384, hidden_dims=(1024, 512), latent_dim=LATENT_DIM):
        super().__init__(input_dim, hidden_dims, latent_dim)


class MetabolomicsEncoder(_MLPEncoder):
    """225 -> 256 -> 2x128"""

    def __init__(self, input_dim=225, hidden_dims=(256,), latent_dim=LATENT_DIM):
        super().__init__(input_dim, hidden_dims, latent_dim)


def product_of_experts(mu_list, logvar_list, include_prior=True):
    mus = torch.stack(mu_list, dim=0)
    logvars = torch.stack(logvar_list, dim=0)
    precisions = torch.exp(-logvars)
    tau_sum = precisions.sum(dim=0)
    mu_weighted = (mus * precisions).sum(dim=0)
    if include_prior:
        tau_sum = tau_sum + 1.0
    var_joint = 1.0 / tau_sum
    mu_joint = var_joint * mu_weighted
    logvar_joint = torch.log(var_joint)
    return (mu_joint, logvar_joint)


class DualEncoderPoE(nn.Module):

    def __init__(
        self,
        trans_dim=17384,
        meta_dim=225,
        trans_hidden=(1024, 512),
        meta_hidden=(256,),
        latent_dim=LATENT_DIM,
    ) -> None:
        super().__init__()
        self.trans_encoder = TranscriptomicsEncoder(trans_dim, trans_hidden, latent_dim)
        self.meta_encoder = MetabolomicsEncoder(meta_dim, meta_hidden, latent_dim)
        self.latent_dim = latent_dim

    def forward(self, x_trans, x_meta):
        mu_t, logvar_t = self.trans_encoder(x_trans)
        mu_m, logvar_m = self.meta_encoder(x_meta)
        mu_joint, logvar_joint = product_of_experts(
            [mu_t, mu_m], [logvar_t, logvar_m], include_prior=True
        )
        return {
            "mu_joint": mu_joint,
            "logvar_joint": logvar_joint,
            "mu_trans": mu_t,
            "logvar_trans": logvar_t,
            "mu_meta": mu_m,
            "logvar_meta": logvar_m,
        }


if __name__ == "__main__":
    torch.manual_seed(0)
    B = 32

    model = DualEncoderPoE()
    x_trans = torch.randn(B, 17384)
    x_meta = torch.randn(B, 225)
    out = model(x_trans, x_meta)

    for k, v in out.items():
        print(f"{k}: {tuple(v.shape)}")

    var_joint = out["logvar_joint"].exp()
    var_t = out["logvar_trans"].exp()
    var_m = out["logvar_meta"].exp()

    n_params = sum(p.numel() for p in model.parameters())
    print(f"\nTotal params: {n_params:,}")

    passed = (
        out["mu_joint"].shape == (B, 128)
        and out["logvar_joint"].shape == (B, 128)
        and torch.all(var_joint <= var_t + 1e-5)
        and torch.all(var_joint <= var_m + 1e-5)
        and torch.isfinite(out["mu_joint"]).all()
        and torch.isfinite(out["logvar_joint"]).all()
    )
    print(f"\n{'PASS' if passed else 'FAIL — check output above.'}")
