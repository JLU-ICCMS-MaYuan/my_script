#!/usr/bin/env python3

import os

for root, dirs, files in os.walk("."):
    if "tmp" in dirs:
        os.system("rm -rf tmp")