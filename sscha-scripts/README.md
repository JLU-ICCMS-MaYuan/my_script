## 流程大纲

### dyn文件改名: 222q网格计算出来的4个dyn文件分别改名为
mv CeSc2H24.dyn1 dyn_pop0_1
mv CeSc2H24.dyn2 dyn_pop0_2
mv CeSc2H24.dyn3 dyn_pop0_3
mv CeSc2H24.dyn4 dyn_pop0_4
dyn_pop0_1 代表第0代的第1个动力学矩阵，第0代也就是初始代

### 1_CreatEns.py生成随机结构
1_CreatEns.py是生成自洽文件的模板，具体到不同的体系，8_ReStart.sh会相应的生成R1_CreatEns.py, 通过运行R1_CreatEns.py生成随机结构。

R1_CreatEns.py将生成的随机结构存放在data_ensemble_manual中，分别是scf_population1_1.dat, ..., scf_population1_100.dat。
scf_population1_1.dat 代表第1代随机结构的第1个随机结构。总共有100个随机结构，随机结构的总数是由8_ReStart.sh中ENS_NUM控制的。

### 2_CreatQEInput.py生成qe自洽输入文件
2_CreatQEInput.py是生成qe自洽输入文件的模板，具体到不同的体系，8_ReStart.sh会相应的生成R2_CreatQEInput.py，通过运行R2_CreatQEInput.py生成全部随机结构qe自洽输入文件。

R2_CreatQEInput.py获取data_ensemble_manual中全部的随结构scf_population1_*.dat, 然后在run_calculation中放置全部随机结构qe自洽输入文件。espresso_run_1.pwi, ..., espresso_run_100.pwi

注意事项：espresso_run_*.pwi 是qe的自洽输入文件，其中的 kpoints 需要测试，最大误差为0.1Ry，每个自洽用时最长为20min， 即：不需要保证k网格是收敛的参数，只要保证在0.1Ry能量误差的精度范围内用时尽可能少即可。kpoints 参数的修改不需要手动进各个 espresso_run_*.pwi 文件修改，只需要在 2_CreatQEInput.py 中修改，它会自动同步生成 R2_CreatQEInput.py. 

### 3_Creat_Sub.py生成提交任务的脚本
这里是需要使用者修改参数的部分，针对你使用的机器和环境，进行修改。具体就是：header变量的内容和29行qe运行命令

如果需要计算总共100个qe的自洽，您想用20个任务，每个任务5个自洽完成，那么您直接更改3_Creat_Sub.py中的这两个值就可以。
```python
n_sub = 20
n_pw=5
```
他会帮您自动生成提交任务的脚本

### 4_SubAllJobs 提交运算 全部随机结构的自洽.

这里是需要使用者修改参数的部分, 具体来说：修改相应作业脚本的提交命令

### 5_CheckFinish.py检查计算情况

在等待计算的过程中可以用5_CheckFinish.py检查计算情况，有多少算完的，有多少还没算的。如果全部交上去的计算任务算完，但是有部分报错或出问题，可以用5B_SubUnfinish.py重新提交。

