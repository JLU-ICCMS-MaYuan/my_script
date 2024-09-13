#!/bin/bash

/lustre/home/h240012/soft/deepmd-kit/bin/python3.10 calypso_run_opt.py || /lustre/home/h240012/soft/deepmd-kit/bin/python3.10 calypso_check_outcar.py > log 2>&1
