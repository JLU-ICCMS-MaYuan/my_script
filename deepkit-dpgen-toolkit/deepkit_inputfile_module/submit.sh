#!/bin/bash

/work/home/may/miniconda3/envs/cage/bin/python calypso_run_opt.py || /work/home/may/miniconda3/envs/cage/bin/python calypso_check_outcar.py > log 2>&1
