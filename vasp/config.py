from argparse import ArgumentParser

class config:
    def __init__(self, args:ArgumentParser):
        self.args = args
    def read_config(self):
        config = {}
        config["input_file_path"] = self.args.input_file_path

        config["press"]           = self.args.press
        config["presses"]         = self.args.presses
        
        config["work_path"]       = self.args.work_path
        config["submit_job_system"] = self.args.submit_job_system

        config["vasp_workflow_type"] = self.args.vasp_workflow
        for other_arg in self.args.more_args:
            # arg_name, value = other_arg.split("=")
            arg_name, value = other_arg.split("=", 1) # 代表只在第一个等号出现的位置劈裂字符串
            config[arg_name] = value
    
        return config