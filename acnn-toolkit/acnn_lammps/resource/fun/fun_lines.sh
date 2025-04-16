#!/bin/bash

echo_line1() {
    local line=$(printf "%-80s" "-")
    echo "${line// /-}"
}
