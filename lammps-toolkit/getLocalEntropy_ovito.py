#!/usr/bin/env python3
import os
import argparse
import numpy as np

import matplotlib.pyplot as plt

from ovito.data import CutoffNeighborFinder, DataCollection
from ovito.modifiers import ComputePropertyModifier
from ovito.io import import_file, export_file


def modify(frame: int, data: DataCollection, cutoff = 5.0, sigma = 0.2, use_local_density = True, compute_average = True, average_cutoff = 5.0):
    # Validate input parameters:
    assert(cutoff > 0.0)
    assert(sigma > 0.0 and sigma < cutoff)
    assert(average_cutoff > 0)

    # Show message in OVITO's status bar:
    yield 'Calculating local entropy'

    # Overall particle density:
    global_rho = data.particles.count / data.cell.volume

    # Initialize neighbor finder:
    finder = CutoffNeighborFinder(cutoff, data)

    # Create output array for local entropy values
    local_entropy = np.empty(data.particles.count)

    # Number of bins used for integration:
    nbins = int(cutoff / sigma) + 1

    # Table of r values at which the integrand will be computed:
    r = np.linspace(0.0, cutoff, num=nbins)
    rsq = r**2

    # Precompute normalization factor of g_m(r) function:
    prefactor = rsq * (4 * np.pi * global_rho * np.sqrt(2 * np.pi * sigma**2))
    prefactor[0] = prefactor[1] # Avoid division by zero at r=0.

    # Iterate over input particles:
    for particle_index in range(data.particles.count):
        yield particle_index / data.particles.count

        # Get distances r_ij of neighbors within the cutoff range.
        r_ij = finder.neighbor_distances(particle_index)

        # Compute differences (r - r_ji) for all {r} and all {r_ij} as a matrix.
        r_diff = np.expand_dims(r, 0) - np.expand_dims(r_ij, 1)

        # Compute g_m(r):
        g_m = np.sum(np.exp(-r_diff**2 / (2.0*sigma**2)), axis=0) / prefactor

        # Estimate local atomic density by counting the number of neighbors within the
        # spherical cutoff region:
        if use_local_density:
            local_volume = 4/3 * np.pi * cutoff**3
            rho = len(r_ij) / local_volume
            g_m *= global_rho / rho
        else:
            rho = global_rho

        # Compute integrand:
        integrand = np.where(g_m >= 1e-10, (g_m * np.log(g_m) - g_m + 1.0) * rsq, rsq)

        # Integrate from 0 to cutoff distance:
        local_entropy[particle_index] = -2.0 * np.pi * rho * np.trapz(integrand, r)

    # Output the computed per-particle entropy values to the data pipeline.
    data.particles_.create_property('Entropy', data=local_entropy)

    # If requested by the user, perform the spatial averaging of the local entropy value using
    # the built-in ComputePropertyModifier. The math expression computes the sum of all entropy values
    # of the particles within the local neighborhood, divided by the number of particles in the neighborhood.
    if compute_average:
        data.apply(ComputePropertyModifier(
            output_property = 'Entropy',
            operate_on = 'particles',
            cutoff_radius = average_cutoff,
            expressions = ['Entropy / (NumNeighbors + 1)'],
            neighbor_expressions = ['Entropy / (NumNeighbors + 1)']))
        

# pipeline = import_file("XDATCAR")
parser = argparse.ArgumentParser(description="Extract a specific frame from XDATCAR and save as POSCAR.")

parser.add_argument("-i", "--filename", default="XDATCAR", help="Path to the XDATCAR file (default: 'XDATCAR').")
args = parser.parse_args()

pipeline = import_file(args.filename)
pipeline.modifiers.append(modify)
pipeline.compute()
# print(pipeline.compute().particles) # output content {'Position': Property('Position'), 'Particle Type': Property('Particle Type'), 'Entropy': Property('Entropy')}
# print(list(pipeline.compute().particles['Entropy']))


# 假设 'Entropy' 数据已经提取为列表 entropy_data
entropy_data = list(pipeline.compute().particles['Entropy'])

# 绘制直方图
plt.figure(figsize=(8, 6))
plt.hist(entropy_data, bins=100, density=True, color='skyblue', edgecolor='black', alpha=0.7)

# 设置标题和标签
plt.title('Entropy Distribution', fontsize=16)
plt.xlabel('Entropy', fontsize=14)
plt.ylabel('Frequency', fontsize=14)

# 网格和显示
plt.grid(axis='y', linestyle='--', alpha=0.7)
# plt.tight_layout()
plt.savefig(os.path.basename(args.filename)+'.png')
# for idx, e in enumerate(pipeline.compute().particles['Entropy']):
#     print(idx+1, e)
# export_file(pipeline, "entropy.dat", "txt/attr", key="Entropy")
# print(pipeline.compute().particles['Entropy / (NumNeighbors + 1)'])