#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import re

def parse_filband(feig, npl=10):
    f = open(feig, 'r')
    lines = f.readlines()

    header = lines[0].strip()
    line = header.strip('\n')
    shape = re.split('[,=/]', line)
    nbnd = int(shape[1])
    nks = int(shape[3])
    eig = np.zeros((nks, nbnd), dtype=np.float32)

    dividend = nbnd
    divisor = npl
    div = nbnd // npl + 1 if nbnd % npl == 0 else nbnd // npl + 2
    kinfo = []
    for index, value in enumerate(lines[1:]):
        value = value.strip(' \n')
        quotient = index // div
        remainder = index % div

        if remainder == 0:
            kinfo.append(value)
        else:
            value = re.split('[ ]+', value)
            a = (remainder - 1) * npl
            b = a + len(value)
            eig[quotient][a:b] = value

    f.close()
    return eig, nbnd, nks, kinfo

def match_pattern_with_line_numbers(pattern, text):
    line_number = 1
    for line in text.splitlines():
        match = re.search(pattern, line)
        if match:
            yield line_number, match
        line_number += 1

def get_sum_orbit(orb, ntype, ielem, i_soc):
    iorb = np.zeros([ntype, ], dtype=np.int32)  # number of projectors for each element
    for i in range(ntype):
        iorb[i] = len(orb[i])
    lorb = np.zeros([ntype, ], dtype=np.int32)  # number of local orbital for each element
    for i in range(ntype):
        for j in orb[i]:
            if j == 's':
                lorb[i] += 1 * i_soc
            elif j == 'p':
                lorb[i] += 3 * i_soc
            elif j == 'd':
                lorb[i] += 5 * i_soc
            elif j == 'f':
                lorb[i] += 7 * i_soc
            else:
                print("unexpected: ", j)
                assert False

    return np.dot(ielem, lorb)

def extract_k_points_and_labels(kpoints_file):
    """
    Extracts the k-points and corresponding labels from the 'eleband.in' file.
    Ensures that the first x_tick is 0 and computes the number of points between high-symmetry points.
    """
    with open(kpoints_file, 'r') as f:
        lines = f.readlines()

    k_points = []
    labels = []
    found_kpoints_section = False
    num_kpoints_between = []  # To store the number of points between each pair of high-symmetry points

    for line in lines:
        if 'K_POINTS {crystal_b}' in line:
            found_kpoints_section = True
        elif found_kpoints_section:
            if line.strip().isdigit():
                continue  # Skip the number of k-points
            elif "!" in line:
                parts = line.split("!")
                label = parts[1].strip()
                num_points = int(parts[0].split()[-1])  # The number of points between high-symmetry points
                num_kpoints_between.append(num_points)
                labels.append(label)
    else:
        num_kpoints_between.pop()
        num_kpoints_between[-1] = num_points + 1


    # Start x_ticks at 0
    x_ticks = [0]
    
    # Add cumulative k-point distances based on the number of points between them
    for i, k_dist in enumerate(num_kpoints_between):
        x_ticks.append(x_ticks[-1] + num_kpoints_between[i])
    print(x_ticks)
    print(labels)

    return x_ticks, labels


def draw_proj_band(proj_file, bd_file, fig_file, ielem, label, color, x_ticks, x_labels, eig_ref, oo):

    eig, nbnd, nks, kinfo = parse_filband(bd_file)

    #eig_ref = max(eig[:, nband - 1])  # VBM for insulators

    ymin = -5  # plot range for y-axis
    ymax = 5
    lw = 0.5  # line width

    scale = 90.0

    ntype = len(ielem)  # number of atom types
    nat = sum(ielem)
    nline_io_header = 0  # line number at '    F    F' or '    T    T'
    with open(proj_file, 'r') as f:
        for line_number in range(nat + ntype + 6 + 3):
            line = f.readline()
            pattern = r"    [TF]    [TF]"
            for i, match in match_pattern_with_line_numbers(pattern, line):
                if match.group()[9] == 'F':
                    i_soc = 1
                elif match.group()[9] == 'T':
                    i_soc = 2
                nline_io_header = line_number + 1
    assert nline_io_header

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    F = plt.gcf()
    F.set_size_inches([6, 4.5])
    p1 = plt.subplot(1, 1, 1)

    # draw bands
    for i in range(nbnd):
        line1 = plt.plot(np.arange(0, nks), eig[:, i] - eig_ref, color='grey', linewidth=lw)

    # draw vertical lines, positions specified by list like [0, 30, 60, ...]
    for vline in x_ticks:
        plt.axvline(x=vline, ymin=ymin, ymax=ymax, linewidth=lw, color='black')

    # Draw horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='--', linewidth=lw)

    with open(proj_file, 'r') as filproj:
        for i in range(nline_io_header - 1):
            l = filproj.readline()  # l : natomwfc, nkstot, nbnd

        nlorb = int(l.split()[0])
        if int(l.split()[1]) != nks or int(l.split()[2]) != nbnd:
            print("Warning: dimension mismatch between proj_file and bd_file!")

        check_nlorb = False
        if check_nlorb:
            orb = [['s', 'p', 'd'], ['s', 'p']]  # projectors for each element
            _nlorb = get_sum_orbit(orb, ntype, ielem, i_soc)
            assert _nlorb == nlorb

        pjmat = np.zeros([nlorb, nks, nbnd], dtype=np.float32)

        _ = filproj.readline()  # skip '    F    F' line
        for i in range(nlorb):
            _ = filproj.readline()  # skip orbit header
            for j in range(nks):
                for k in range(nbnd):
                    pjmat[i, j, k] = float(filproj.readline().split()[2])

    for i in range(len(oo)):  # once plot a type
        plt.scatter(-1, ymin - 1, 20, c=color[i], alpha=0.5, label=label[i], marker='.', edgecolor='none')
        for k in range(nbnd):  # once plot a band
            s_of_o = np.zeros([nks, ], dtype=np.float32)  # size of dots for all kpoints in a band
            for j in oo[i]:
                s_of_o[:] += pjmat[j, :, k]

            plt.scatter(np.arange(0, nks), eig[:, k] - eig_ref, s=scale * s_of_o, c=color[i], alpha=0.5, marker='.',
                        edgecolor='none')

    plt.xlim([0, nks - 1])
    plt.ylim([ymin, ymax])
    plt.ylabel(r'E (eV)', fontsize=16)
    plt.xticks(x_ticks, x_labels)

    plt.subplots_adjust(left=0.20, right=0.75, top=0.95, bottom=0.1)
    p1.legend(scatterpoints=1, numpoints=1, markerscale=2.0, bbox_to_anchor=(1.05, 1), loc='upper left',
              borderaxespad=0.)

    plt.savefig(fig_file, dpi=500)

if __name__ == '__main__':
    info = '''
Some notes:
I'll give you an example of using this scripts:
python plot_band_for_project.py pjband.png --ielem 4 14 --label 'Nb d' 'H s' --color r g --efermi 23.4483 --oo '4 5 6 7 8 17 18 19 20 21 30 31 32 33 34 43 44 45 46 47' '52 53 54 55 56 57 58 59 60 61 62 63 64 65'

The most important part is to specify the number of orbis for each element.

You can get it with the following code:
grep '[a-zA-Z]'  elefatbandpro.projwfc_up 

grep '[a-zA-Z]'  elefatbandpro.projwfc_up | grep 1S | awk '{print $1-1}' | xargs
52 53 54 55 56 57 58 59 60 61 62 63 64 65 
There are 14 H's in the system, and each H has 1 1s orbital, so there are 14 numbers

grep '[a-zA-Z]'  elefatbandpro.projwfc_up | grep "Nb  4D" | awk '{print $1-1}' | xargs
4 5 6 7 8 17 18 19 20 21 30 31 32 33 34 43 44 45 46 47
There are 4 Nb's in the system, and each Nb's 4d orbital has 5 magnetic quantum number components, so there are 20 numbers
'''
    print(info)
    parser = argparse.ArgumentParser(description="Draw projected band structure")
    parser.add_argument("fig_file", type=str, help="Path to the output figure file")
    parser.add_argument("--ielem", type=int, nargs='+', required=True, help="Number of atoms for each element (e.g., 4 14)")
    parser.add_argument("--label", type=str, nargs='+', required=True, help="Element and orbit labels (e.g., 'Nb s' 'Nb p' 'Nb d' 'H s')")
    parser.add_argument("--color", type=str, nargs='+', required=True, help="Colors for each orbit (e.g., 'r' 'g' 'b')")
    parser.add_argument("--efermi", type=float, required=True, help="the fermi energy level")
    parser.add_argument("--oo", type=str, nargs='+', required=True, help="Orbital indices wrapped in quotes, e.g., '9 22 35 48'")

    args = parser.parse_args()

    # Convert oo (string representation) to a list of lists of integers
    oo = [list(map(int, group.split())) for group in args.oo]

    # Default file paths
    proj_file = 'elefatbandpro.projwfc_up'
    bd_file = 'elebanddata.dat'
    kpoints_file = 'eleband.in'

    # Extract x_ticks and x_labels from eleband.in
    x_ticks, x_labels = extract_k_points_and_labels(kpoints_file)


    draw_proj_band(proj_file, bd_file, args.fig_file, args.ielem, args.label, args.color, x_ticks, x_labels, args.efermi, oo)
    print("done")
