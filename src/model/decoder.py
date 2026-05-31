import torch
import torch.nn as nn

LATENT_DIM = 128


class _MLPDecoder(nn.Module):

    def __init__(self, latent_dim, hidden_dims, output_dim):
        super().__init__()
        layers = []
        prev = latent_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.LayerNorm(h), nn.GELU()]
            prev = h
        self.backbone = nn.Sequential(*layers)
        self.head = nn.Linear(prev, output_dim)
        self.soft_assignment = None
        self.output_dim = output_dim

    def forward(self, z):
        if self.soft_assignment is not None:
            raise NotImplementedError("KEGG soft assignment activated later.")
        return self.head(self.backbone(z))


class TranscriptomicsDecoder(_MLPDecoder):
    """128 → 512 → 1024 → 17,384"""

    def __init__(
        self, latent_dim=LATENT_DIM, hidden_dims=(512, 1024), output_dim=17384
    ):
        super().__init__(latent_dim, hidden_dims, output_dim)


class MetabolomicsDecoder(_MLPDecoder):
    """128 → 256 → 225"""

    def __init__(self, latent_dim=LATENT_DIM, hidden_dims=(256,), output_dim=225):
        super().__init__(latent_dim, hidden_dims, output_dim)


if __name__ == "__main__":
    torch.manual_seed(0)
    B = 32

    t_decoder = TranscriptomicsDecoder()
    m_decoder = MetabolomicsDecoder()

    z = torch.randn(B, 128)
    x_trans_hat = t_decoder(z)
    x_meta_hat = m_decoder(z)
    n_params = sum(p.numel() for p in t_decoder.parameters())
    print(f"\nTotal params: {n_params:,}")
    passed = (
        x_trans_hat.shape == (B, 17384)
        and x_meta_hat.shape == (B, 225)
        and torch.isfinite(x_trans_hat).all()
        and torch.isfinite(x_meta_hat).all()
        and t_decoder.soft_assignment is None
        and m_decoder.soft_assignment is None
    )
    print(f"\n{'PASS' if passed else 'FAIL — check output above.'}")
