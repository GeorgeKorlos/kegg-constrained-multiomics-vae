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


def compute_eta(loss_trans, loss_meta, shared_params, eps=1e-8):
    g_trans = torch.autograd.grad(
        loss_trans, shared_params, retain_graph=True, allow_unused=True
    )
    g_meta = torch.autograd.grad(
        loss_meta, shared_params, retain_graph=True, allow_unused=True
    )

    def grad_norm(grads):
        sq = [(g**2).sum() for g in grads if g is not None]
        if not sq:
            return torch.zeros((), device=loss_trans.device)
        return torch.sqrt(torch.stack(sq).sum())

    norm_trans = grad_norm(g_trans)
    norm_meta = grad_norm(g_meta)
    return (norm_trans / (norm_meta + eps)).detach()
