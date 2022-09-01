import os
import logging
from argparse import ArgumentParser
from pathlib import Path
from pprint import pprint


from config import config
from specify_wyckoffs import specify_wyckoffs
from substitution import substitution
from pso import pso

class generator_methods:

    def __init__(self, args: ArgumentParser) -> None:
        
        self.config_d = config(args).config_d

        if self.config_d["mode"] == "specifywps":
            specify_wyckoffs.init_from_config(self.config_d)

        if self.config_d["mode"] == "substitution":
            substitution.init_from_config(self.config_d)

        if self.config_d["mode"] == "pso":
            pso.init_from_config(self.config_d)
