import torch
import torch.nn as nn

LATENT_DIM = 128

from src.model.encoder import DualEncoderPoE
from src.model.decoder import TranscriptomicsDecoder, MetabolomicsDecoder


def reparameterize(mu, logvar, sample=True):
    if not sample:
        return mu
    std = torch.exp(0.5 * logvar)
    eps = torch.randn_like(std)
    return mu + eps * std


class MultiOmicsVAE(nn.Module):

    def __init__(
        self,
        trans_dim=17384,
        meta_dim=225,
        latent_dim=LATENT_DIM,
        trans_hidden=(1024, 512),
        meta_hidden=(256,),
    ):
        super().__init__()
        self.encoder = DualEncoderPoE(
            trans_dim, meta_dim, trans_hidden, meta_hidden, latent_dim
        )
        self.trans_decoder = TranscriptomicsDecoder(
            latent_dim, tuple(reversed(trans_hidden)), trans_dim
        )
        self.meta_decoder = MetabolomicsDecoder(
            latent_dim, tuple(reversed(meta_hidden)), meta_dim
        )

    def forward(self, x_trans, x_meta, sample=True):
        enc = self.encoder(x_trans, x_meta)
        z = reparameterize(enc["mu_joint"], enc["logvar_joint"], sample)
        x_trans_hat = self.trans_decoder(z)
        x_meta_hat = self.meta_decoder(z)
        return {**enc, "z": z, "x_trans_hat": x_trans_hat, "x_meta_hat": x_meta_hat}


if __name__ == "__main__":
    torch.manual_seed(0)
    B = 32

    model = MultiOmicsVAE()
    x_trans = torch.randn(B, 17384)
    x_meta = torch.randn(B, 225)
    out = model(x_trans, x_meta, sample=True)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"\nTotal params: {n_params:,}")
    passed = (
        out["x_trans_hat"].shape == (B, 17384)
        and out["x_meta_hat"].shape == (B, 225)
        and out["z"].shape == (B, 128)
        and out["mu_joint"].shape == (B, 128)
    )

    out1 = model(x_trans, x_meta, sample=False)
    out2 = model(x_trans, x_meta, sample=False)
    assert torch.equal(out1["z"], out2["z"])
    assert torch.equal(out1["x_trans_hat"], out2["x_trans_hat"])
    assert torch.equal(out1["x_meta_hat"], out2["x_meta_hat"])
    print(f"\n{'PASS' if passed else 'FAIL — check output above.'}")
