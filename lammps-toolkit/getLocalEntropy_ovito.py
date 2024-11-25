#!/usr/bin/env python3
##### Calculate local entropy #####
#
# This modifier function computes the local pair entropy fingerprint for each particle in the system.
# [See the documentation](manual:modifiers.calculate_local_entropy).

# This is a copy of the file "E:\OVITO Basic\plugins\python\ovito\_extensions\scripts\modifiers\Calculate local entropy.py".
# Feel free to modify the Python code below as needed. Changes will NOT affect the original master file.
import argparse

from ovito.io import import_file, export_file
from ovito.data import CutoffNeighborFinder, DataCollection
from ovito.modifiers import ComputePropertyModifier
from ovito.pipeline import ModifierInterface
import numpy as np

class MyModifier(ModifierInterface):

    cutoff=5.0, 
    sigma=0.2,  
    use_local_density=True, 
    compute_average=True, 
    average_cutoff=5.0

    def modify(self, frame: int, data: DataCollection, cutoff = 5.0, sigma = 0.2, use_local_density = False, compute_average = False, average_cutoff = 5.0):
    # def modify(self, cutoff = 5.0, sigma = 0.2, use_local_density = False, compute_average = False, average_cutoff = 5.0):
        # Validate input parameters:
        # assert(cutoff > 0.0)
        # assert(sigma > 0.0 and sigma < cutoff)
        # assert(average_cutoff > 0)
        print(cutoff, sigma, average_cutoff)
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

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Calculate local entropy for particles.")
    
    # Add arguments
    parser.add_argument('-f', '--frame', type=int, help="Frame number")
    parser.add_argument('-d', '--data', type=str, help="Path to the input data file")
    parser.add_argument('-r', '--cutoff', type=float, default=5.0, help="Cutoff distance")
    parser.add_argument('-s', '--sigma', type=float, default=0.2, help="Sigma value for Gaussian function")
    parser.add_argument('-l', '--local_density', action='store_true', help="Use local density for entropy calculation")
    parser.add_argument('-c', '--compute_average', action='store_true', help="Compute spatial average of entropy")
    parser.add_argument('-a', '--average_cutoff', type=float, default=5.0, help="Cutoff distance for averaging")
    parser.add_argument("-o", "--outputname", default='entropy_output.txt', type=str, help="Output filename for entropy values")

    # Parse arguments
    args = parser.parse_args()

    # Import the file (use correct path)
    pipeline = import_file(args.data)
    print("Number of MD frames:", pipeline.num_frames)

    # Add the custom local entropy modifier to the pipeline
    pipeline.modifiers.append(MyModifier(cutoff=args.cutoff, sigma=args.sigma,  use_local_density=args.local_density, compute_average=args.compute_average, average_cutoff=args.average_cutoff))
    # (frame=1, data=pipeline.compute(), cutoff=args.cutoff, sigma=args.sigma,  use_local_density=args.local_density, compute_average=args.compute_average, average_cutoff=args.average_cutoff)
    data = pipeline.compute()
    print(data.particles['Entropy'])
    # export_file(pipeline,args.outputname, "txt/table", key="Entropy")