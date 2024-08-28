#!/bin/bash

for i in `cat log`; do 
    mkdir retry_$i
    cp $i.structures retry_$i -fr
    cp calypso_run_model_devi.py retry_$i
    cd retry_$i
    python calypso_run_model_devi.py --all-models ../../gen_stru_analy.000/graph.000.pb  ../../gen_stru_analy.000/graph.001.pb  ../../gen_stru_analy.000/graph.002.pb  ../../gen_stru_analy.000/graph.003.pb  --type_map La Sc H &
    cd ..
done
