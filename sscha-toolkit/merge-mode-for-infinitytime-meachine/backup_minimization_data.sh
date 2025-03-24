#!/bin/bash

begin_idx=$1
end_idx=$2

cp minimization_data.dat    minimization_data_${begin_idx}_${end_idx}.dat
cp minimization_data.freqs  minimization_data_${begin_idx}_${end_idx}.freqs
cp sscha.out                sscha_${begin_idx}_${end_idx}.out
cp save_filename.dat        save_filename_${begin_idx}_${end_idx}.dat
