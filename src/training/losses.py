import torch

F_TRANS = 17384
F_META = 225
LATENT_DIM = 128


def reconstruction_loss_trans(x_hat, x):
    return ((x_hat - x) ** 2).sum(dim=-1).mean() / F_TRANS


def reconstruction_loss_meta(x_hat, x):
    return ((x_hat - x) ** 2).sum(dim=-1).mean() / F_META


def kl_loss(mu, logvar):
    kl = 0.5 * (mu**2 + logvar.exp() - logvar - 1.0)
    return kl.sum(dim=-1).mean() / LATENT_DIM


def compute_eta(loss_trans, loss_meta, z, eps=1e-8):
    g_trans = torch.autograd.grad(loss_trans, z, retain_graph=True)[0]
    g_meta = torch.autograd.grad(loss_meta, z, retain_graph=True)[0]
    eta = g_trans.norm() / (g_meta.norm() + eps)
    return eta.detach()
