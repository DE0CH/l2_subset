# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
from code_snippet_FC import mixture_of_gaussians, sample_from_mixture
import torch
import numpy as np
import matplotlib.pyplot as plt
import re
import argparse

# %%
def plot_contour_with_points(points, target_distribution, xlim=(-4, 4), ylim=(-4, 4), levels=50):

    x = np.linspace(xlim[0], xlim[1], 100)
    y = np.linspace(ylim[0], ylim[1], 100)
    X, Y = np.meshgrid(x, y)
    XY = np.stack([X, Y], axis=-1)
    XY_tensor = torch.tensor(XY, dtype=torch.float32)

    log_p, _ = target_distribution(XY_tensor)
    log_p_np = log_p.cpu().detach().numpy()

   
    fig, ax = plt.subplots(figsize=(8, 6))
    contour = ax.contourf(X, Y, log_p_np, levels=levels, cmap="viridis")
    fig.colorbar(contour, ax=ax, label="Log Density")

    # Plot points if provided
    if isinstance(points, torch.Tensor):
        points = points.cpu().detach().numpy()
    ax.scatter(points[:, 0], points[:, 1], color="red", edgecolor="black", label="Sampled Points")

    ax.set_title("Subset Selected Points")
    ax.legend()

    return fig


# %%
def find_last_active_points_line(filename):
    last_active_line = None
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("Active points: "):
                last_active_line = line.strip()  # Remove newline characters
    return last_active_line


# %%
def get_numbers_from_line(line):
    match = re.search(r'Active points:\s*(.*)', line)
    if match:
        numbers_str = match.group(1)  # Get the captured group
        numbers = list(map(int, numbers_str.split()))  # Convert space-separated numbers to integers
        return numbers
    raise ValueError('The line is malformed')

def extract_n_d_from_file(filename):
    pattern = re.compile(r'Running a new restart with seed \d+, n = (\d+), m = \d+, d = (\d+), p = \d+')
    
    with open(filename, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                return int(match.group(1)), int(match.group(2)) # Return the extracted 'n' and 'd' value as an integer
    raise ValueError('The file does not contain the expected pattern')

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('subset_log', type=str, help='Log file of subset selection')
    parser.add_argument('output_file', type=str, help='Output file')
    parser.add_argument('--seed', type=int, default=42)

    args = parser.parse_args()
    n, d = extract_n_d_from_file(args.subset_log)
    points = sample_from_mixture(n, d, seed=args.seed)
    indices = torch.tensor(get_numbers_from_line(find_last_active_points_line(args.subset_log)))
    fig = plot_contour_with_points(points[indices], mixture_of_gaussians)
    fig.savefig(args.output_file)

# %%
if __name__ == '__main__':
    main()
