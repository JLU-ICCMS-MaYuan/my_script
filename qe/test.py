from pathlib import Path
import os
files_underdir = os.listdir("/home/mayuan/mycode/my_script/test/test_qe/200.0")
for file in files_underdir:
    print(Path(file).suffix); print(type(Path(file).suffix))
    if Path(file).suffix == ".dyn0":
        print("----------------------------------------------")
else:
    print(f"There is no suffix named dyn0") 