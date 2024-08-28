import os  
import glob  

pattern = 'task.*.000'  
tasks = glob.glob(pattern)
tasks.sort(key=lambda x: int(x.split('.')[1])) 
# 遍历文件并处理  
idx = 0
with open('Model_Devi.out', 'w') as f1:
    f1.write("#\n#       step         max_devi_v         min_devi_v         avg_devi_v         max_devi_f         min_devi_f         avg_devi_f             devi_e            min_dis\n")
    for task in tasks:  
        # 这里只打印文件名作为示例  
        # 你可以添加更多的代码来读取文件内容并进行处理  
        model_devi_out = os.path.join(task, 'model_devi.out')
        if os.path.exists(model_devi_out):
            # 示例：读取文件内容  
            with open(model_devi_out, 'r') as f2:  
                lines = f2.readlines()[2:]
            for line in lines:
                line = line.split()
                if len(line)==9:
                    line[0] = str(idx)
                    idx += 1
                    f1.write('       '.join(line)+'\n')
                else:
                    print(model_devi_out)
                    break
        else:
            print(model_devi_out)