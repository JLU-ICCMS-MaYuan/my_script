#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert crystal matrix to 6 lattice parameters (a, b, c, α, β, γ)
Supported input formats: POSCAR/CONTCAR (VASP), or direct 3x3 matrix input
"""

import numpy as np
import sys
import math

def read_lattice_from_poscar(poscar_path):
    """Read lattice vectors from VASP POSCAR/CONTCAR"""
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    
    # Skip comment line and scaling factor
    scale = float(lines[1].strip())
    lattice_vectors = []
    for line in lines[2:5]:
        lattice_vectors.append([float(x) for x in line.split()[:3]])
    
    return np.array(lattice_vectors) * scale

def matrix_to_parameters(lattice_matrix):
    """
    Convert 3x3 lattice matrix to 6 parameters (a, b, c, α, β, γ)
    α (alpha): angle between b and c (in degrees)
    β (beta): angle between a and c
    γ (gamma): angle between a and b
    """
    a = np.linalg.norm(lattice_matrix[0])
    b = np.linalg.norm(lattice_matrix[1])
    c = np.linalg.norm(lattice_matrix[2])
    
    alpha = np.degrees(np.arccos(np.dot(lattice_matrix[1], lattice_matrix[2]) / (b * c)))
    beta  = np.degrees(np.arccos(np.dot(lattice_matrix[0], lattice_matrix[2]) / (a * c)))
    gamma = np.degrees(np.arccos(np.dot(lattice_matrix[0], lattice_matrix[1]) / (a * b)))
    
    return a, b, c, alpha, beta, gamma

def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_crystalmatrix2crystalparameter.py [POSCAR/CONTCAR] or")
        print("       python convert_crystalmatrix2crystalparameter.py a11 a12 a13 a21 a22 a23 a31 a32 a33")
        sys.exit(1)
    
    try:
        if len(sys.argv) == 2:  # POSCAR path provided
            poscar_path = sys.argv[1]
            lattice_matrix = read_lattice_from_poscar(poscar_path)
            print(f"\nReading lattice from: {poscar_path}")
        else:  # Direct matrix elements provided
            elements = [float(x) for x in sys.argv[1:]]
            if len(elements) != 9:
                raise ValueError("Need exactly 9 elements for 3x3 matrix")
            lattice_matrix = np.array(elements).reshape(3,3)
            print("\nUsing direct lattice matrix input:")
        
        print("\nLattice matrix (Å):")
        for row in lattice_matrix:
            print(f"  {row[0]:10.6f} {row[1]:10.6f} {row[2]:10.6f}")
        
        a, b, c, alpha, beta, gamma = matrix_to_parameters(lattice_matrix)
        
        print("\nLattice parameters:")
        print(f"  a     = {a:10.6f} Å")
        print(f"  b     = {b:10.6f} Å")
        print(f"  c     = {c:10.6f} Å")
        print(f"  α     = {alpha:10.6f}°")
        print(f"  β     = {beta:10.6f}°")
        print(f"  γ     = {gamma:10.6f}°")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()