#!/bin/bash

for it in RELAX/IT*; do
    [ -d "$it" ] || continue
    cnt=$(find "$it" -maxdepth 1 -type f -name 'group_*' -print0 | xargs -0 cat 2>/dev/null | wc -l)
    printf '%s %d\n' "$it" "$cnt"
done

find RELAX -type f -name 'group_*' -print0 | xargs -0 cat | wc -l
