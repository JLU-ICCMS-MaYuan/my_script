#!/bin/bash

file=$1

find ./ -name "${file}.save" | xargs rm -fr
find ./ -type f \( -name "*wfc*" -o -name "*dwf*"  -o -name "*prd*" -o -name "*bar*"  -o -name "*recover*"  -o -name "*mixd*"  \) | xargs rm -rf