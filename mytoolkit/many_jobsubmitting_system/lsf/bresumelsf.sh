#!/bin/bash

for tkid in `awk '{print $2}' tasks`; do
    bresume $tkid
done

