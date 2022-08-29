from pathlib import Path
from configparser import ConfigParser
from argparse import ArgumentParser

class config:
    def __init__(self, args: ArgumentParser) -> None:
        self.args = args
        # init the self.config by using the method self.read_config()
        self.config_d = {}
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

        config_base = config_inputini._sections["base"]
        for key, value in config_base.items():
            try:
                self.config_d[key] = eval(value)
            except (NameError, SyntaxError):
                pass

        
