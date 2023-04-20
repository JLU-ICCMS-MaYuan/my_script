import sys
from pathlib import Path
from configparser import ConfigParser
from argparse import ArgumentParser


class config:
    def __init__(self, args: ArgumentParser) -> None:
        self.args = args
        # init the self.config by using the method self.read_config()
        self.config_d = {}
        self.sections = []
        self.read_args()
        self.read_ini()
 

    def read_args(self):
        self.config_d["work_path"]       = self.args.work_path
        self.config_d["input_file_path"] = Path(self.args.input_file_path)
        self.config_d["generator"]       = self.args.generator        

        for other_args in self.args.more_args:
            args_name, value = other_args.split("=")
            self.config_d[args_name] = value

 
    def read_ini(self):

        config_inputini      = ConfigParser()
        config_inputini.read(self.config_d["input_file_path"])

        if "specifywps" in config_inputini.sections():
            self.sections.append("specifywps")
            for key, value in config_inputini.items("specifywps"):
                try:
                    self.config_d[key] = eval(value)
                except (NameError, SyntaxError):
                    pass

        if "splitwps" in config_inputini.sections():
            self.sections.append("splitwps")
            for key, value in config_inputini.items("splitwps"):
                try:
                    self.config_d[key] = eval(value)
                except (NameError, SyntaxError):
                    pass
        
        if "substitution" in config_inputini.sections():
            self.sections.append("substitution")
            for key, value in config_inputini.items("substitution"):
                try:
                    self.config_d[key] = eval(value)
                except (NameError, SyntaxError):
                    print("there is no {}".format(key))
                    sys.exit(1)
        
        if "pso" in config_inputini.sections():
            self.sections.append("pso")
            for key, value in config_inputini.items("pso"):
                try:
                    self.config_d[key] = eval(value)
                except (NameError, SyntaxError):
                    print("there is no {}".format(key))
                    sys.exit(1)
        
