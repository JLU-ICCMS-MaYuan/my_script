from argparse import ArgumentParser

class config:
    def __init__(self, args:ArgumentParser):
        self.args = args
    def read_config(self):
        config = {}
        config["input_file_path"] = self.args.input_file_path
        config["work_path"]       = self.args.work_path
        config["workflow_type"]   = self.args.tool_flow
        for other_arg in self.args.more_args:
            arg_name, value = other_arg.split("=")
            config[arg_name] = value 
        return config