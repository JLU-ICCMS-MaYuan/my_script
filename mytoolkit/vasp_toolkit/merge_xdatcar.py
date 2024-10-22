#!/bin/bash
import os
import sys

from ase.io.vasp import write_vasp_xdatcar, read_vasp_xdatcar

print("Please input XDATCARs in order you wish")
filenames = sys.argv[1:]

xdatcars = []
for filename in filenames:
    num_configs = os.popen("grep 'Direct configuration' {} | wc -l".format(filename)).read().strip()
    xdatcar = read_vasp_xdatcar(filename, slice(0, int(num_configs)))
    print('filename:{}\nfound={}\nread={}'.format(filename, num_configs, len(xdatcar)))
    xdatcars.extend(xdatcar)

print('total configuartions: {}'.format(len(xdatcars)))
write_vasp_xdatcar("xdatcar_merged", xdatcars)
