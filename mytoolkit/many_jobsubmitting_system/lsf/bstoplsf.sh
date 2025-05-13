#!/bin/bash

for tkid in `awk '{print $1}' tasks`; do
    bstop $tkid
done

