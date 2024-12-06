#!/bin/bash

resource_path="./resource"
iter_sch_path="$resource_path/iter_scheduling_remote.sh"
max_iter=100

set -e
echo "Scheduling tast date: $(date)"
echo "Scheduling tast pid: $$"

for iter in $(seq 0 $max_iter);do
    echo "iter: ================================================ $iter"
    bash $iter_sch_path $resource_path
done
