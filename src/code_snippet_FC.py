import torch
import numpy as np
from util import save_matrix_to_binary

gen = torch.Generator()
gen.manual_seed(42)

def mixture_of_gaussians(X):
    """
    Compute the log-density and gradient of the log-density for a mixture of two bivariate Gaussians.
    """
    # Define means and covariances for the two Gaussians
    mean1 = torch.tensor([-1.5, 0.0])
    cov1 = torch.tensor([[1.0, 0.0], [0.0, 1.0]])
    cov_inv1 = torch.inverse(cov1)
    det_cov1 = torch.det(cov1)
    
    mean2 = torch.tensor([1.5, 0.0])
    cov2 = torch.tensor([[1.0, 0.0], [0.0, 1.0]])
    cov_inv2 = torch.inverse(cov2)
    det_cov2 = torch.det(cov2)
    
    # Compute the difference from the means
    diff1 = X - mean1  # Shape: [nbatch, nsamples, 2]
    diff2 = X - mean2  # Shape: [nbatch, nsamples, 2]

    # Compute the log densities of each Gaussian (up to a constant)
    log_p1 = -0.5 * torch.einsum('bni,ij,bnj->bn', diff1, cov_inv1, diff1) - 0.5 * torch.log(det_cov1)
    log_p2 = -0.5 * torch.einsum('bni,ij,bnj->bn', diff2, cov_inv2, diff2) - 0.5 * torch.log(det_cov2)
    
    # Compute the gradient of the log-densities
    grad_log_p1 = -torch.einsum('ij,bnj->bni', cov_inv1, diff1)
    grad_log_p2 = -torch.einsum('ij,bnj->bni', cov_inv2, diff2)

    # Use hard assignment based on which Gaussian is more likely
    # Create a mask for which component dominates each point
    mask1 = log_p1 > log_p2
    mask2 = ~mask1

    # Combine the log densities using the mask
    log_p_mixture = torch.where(mask1, log_p1, log_p2)
    grad_log_p_mixture = torch.where(mask1.unsqueeze(-1), grad_log_p1, grad_log_p2)
    
    return log_p_mixture, grad_log_p_mixture


# KSD loss to be used within MPMC torch environment
def KSD_loss(X, nbatch, nsamples, dim):
    X = X.view(nbatch, nsamples, dim)
    log_p, grad_log_p = mixture_of_gaussians(X)

    # Apply temperature scaling to the gradients to emphasize high-density areas
    T = 0.5  # Experiment with values like 0.1 or 0.05
    grad_log_p /= T

    alpha = 0.01
    beta = -0.05
    X_expanded = X.unsqueeze(2)
    X_diff = X_expanded - X_expanded.transpose(1, 2)
    norm_sq = torch.sum(X_diff ** 2, dim=-1)
    K = (alpha + norm_sq) ** beta 

    # Compute gradients of the kernel with respect to x and x'
    grad_K = -2 * beta * X_diff * K.unsqueeze(-1) / (alpha + norm_sq).unsqueeze(-1)

    # Hessian trace term (for divergence)
    hessian_trace = -2 * beta * (dim * (alpha + norm_sq) ** (beta - 1) - 
                                 2 * beta * norm_sq * (alpha + norm_sq) ** (beta - 2))

    grad_log_p_expanded = grad_log_p.unsqueeze(2)
    dot_grad_log_p = torch.einsum('bikd,bjkd->bij', grad_log_p_expanded, grad_log_p_expanded)
    dot_grad_K_x = torch.einsum('bikd,bijd->bij', grad_log_p_expanded, grad_K)
    dot_grad_K_y = torch.einsum('bjkd,bijd->bij', grad_log_p_expanded, grad_K)

    stein_kernel = K * dot_grad_log_p + dot_grad_K_x + dot_grad_K_y + hessian_trace
             
    return stein_kernel


nbatch = 1  # Number of batches
nsamples = 200  # Number of samples per batch
dim = 2  # Dimensionality of each sample
m = 100 # m as pass into L2


X = torch.randn(nbatch, nsamples, dim, generator=gen)

v = KSD_loss(X, nbatch, nsamples, dim)

w = v[0].numpy().astype(np.float64)
print(w)
print(sum(sum(w)))

#save_matrix_to_binary(filename, M, n, m, d) saves the matrix M (as a numpy 2D array) into the file format expected by the l2 code. We selected a m*m matrix from a n*n matrix by minimizing the sum of all entries.
save_matrix_to_binary('t.bin', w, nsamples, m)


