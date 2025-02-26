#!/usr/bin/env python3
import torch
import numpy as np
from utils import save_matrix_to_binary, save_points_to_binary
import argparse
import math

def mixture_of_gaussians(X):
    is_batched = (X.dim() == 3)    # Check if input is batched

    # Define Gaussian mixture components
    means = torch.tensor([[-1.5, 0.0], [1.5, 0.0]])    
    covs = torch.stack([torch.eye(2), torch.eye(2)])
    weights = torch.tensor([0.5, 0.5])    

    if is_batched:
            means = means.unsqueeze(0)    
            X_exp = X.unsqueeze(2) 
    else:
            X_exp = X.unsqueeze(1) 

    # Compute differences
    diff = X_exp - means    

    det_covs = torch.det(covs) 
    inv_covs = torch.inverse(covs)    

    mahalanobis = torch.einsum('...mi,mij,...mj->...m', diff, inv_covs, diff)    

    log_coeff = torch.log(weights) - 0.5 * (2 * torch.log(torch.tensor(2 * math.pi)) + torch.log(det_covs))
    log_exponent = -0.5 * mahalanobis    
    log_p_comp = log_coeff + log_exponent    

    max_log_p = torch.max(log_p_comp, dim=-1, keepdim=True)[0]    
    log_p = max_log_p.squeeze(-1) + torch.log(torch.sum(torch.exp(log_p_comp - max_log_p), dim=-1))    

    # Compute gradient of log density
    weighted_gradients = weights * torch.exp(log_p_comp - log_p.unsqueeze(-1))    
    grad_log_p = -torch.einsum('...m,...mi,mij->...i', weighted_gradients, diff, inv_covs) 

    return log_p, grad_log_p

def KSD_loss_RBF(X, nbatch, nsamples, dim):
    X = X.view(nbatch, nsamples, dim)

    log_p, grad_log_p = mixture_of_gaussians(X)

    # Compute pairwise squared Euclidean distances
    X_expanded = X.unsqueeze(2)    
    X_diff = X_expanded - X_expanded.transpose(1, 2)    
    norm_sq = torch.sum(X_diff ** 2, dim=-1)    

    # Median heuristic for bandwidth selection
    median_dist = torch.median(norm_sq[norm_sq > 0])    
    h = torch.sqrt(0.5 * median_dist / torch.log(torch.tensor(nsamples + 1.0)))

    # RBF kernel computation
    K = torch.exp(-norm_sq / (2 * h**2))

    # Compute the gradient of the kernel
    grad_K = -K.unsqueeze(-1) * X_diff / (h**2)    

    # Compute Hessian trace
    hessian_trace = K * (dim / (h**2) - norm_sq / (h**4))    

    # Compute KSD components
    grad_log_p_expanded = grad_log_p.unsqueeze(2)    
    dot_grad_log_p = torch.einsum('bikd,bjkd->bij', grad_log_p_expanded, grad_log_p_expanded)
    dot_grad_K_x = torch.einsum('bikd,bijd->bij', grad_log_p_expanded, grad_K)
    dot_grad_K_y = torch.einsum('bjkd,bijd->bij', grad_log_p_expanded, grad_K)

    # Compute final Stein kernel
    stein_kernel = K * dot_grad_log_p + dot_grad_K_x + dot_grad_K_y + hessian_trace

    return stein_kernel


def sample_from_mixture(n_samples, dim, seed=None):
    """Generates samples from the target Gaussian mixture model."""
    if seed is not None:
        torch.manual_seed(seed)
    means = torch.tensor([[-1.5, 0.0], [1.5, 0.0]])    # (nMix, nDim)
    covs = torch.stack([torch.eye(2), torch.eye(2)])    # (nMix, nDim, nDim)
    weights = torch.tensor([0.5, 0.5])    # (nMix,)

    component_indices = torch.multinomial(weights, num_samples=n_samples, replacement=True)

    sampled_points = torch.empty((n_samples, dim))
    for i, comp in enumerate(component_indices):
            sampled_points[i] = torch.distributions.MultivariateNormal(means[comp], covs[comp]).sample()

    return sampled_points

if __name__ == '__main__':
    nbatch = 1  # Number of batchess
    parser = argparse.ArgumentParser()
    parser.add_argument('nsamples', type=int)
    parser.add_argument('dim', type=int)
    parser.add_argument('m', type=int)
    parser.add_argument('matrix_file', type=str)
    parser.add_argument('points_file', type=str)
    parser.add_argument('--seed', type=int, default=42)

    args = parser.parse_args()

    nsamples = args.nsamples # Number of samples per batch
    dim = args.dim # Dimensionality of each sample
    m = args.m # m as pass into L2

    X = sample_from_mixture(nsamples, dim, seed=42)

    v = KSD_loss_RBF(X, nbatch, nsamples, dim)

    w = v[0].numpy().astype(np.float64)

    #save_matrix_to_binary(filename, M, n, m) saves the matrix M (as a numpy 2D array) into the file format expected by the l2 code. We selected a m*m matrix from a n*n matrix by minimizing the sum of all entries.
    save_matrix_to_binary(args.matrix_file, w, nsamples, m)
    save_points_to_binary(args.points_file, X.numpy().astype(np.float64))


