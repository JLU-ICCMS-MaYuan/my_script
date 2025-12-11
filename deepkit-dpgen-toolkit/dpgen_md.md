# dpgen+lammps分子动力学

##  <span style="color:red"> 1. deepgen+lammps 的使用注意事项

### <span style="color:yellow"> dpgen 跑了很多代，但是一直不进行FP计算
```shell
2025-11-13 01:29:22,467 - INFO : system 000 candidate :      5 in    620   0.81 %
2025-11-13 01:29:22,468 - INFO : system 000 failed    :    595 in    620  95.97 %
2025-11-13 01:29:22,468 - INFO : system 000 accurate  :     20 in    620   3.23 %
2025-11-13 01:29:22,489 - INFO : system 000 accurate_ratio:   0.0323    thresholds: 1.0000 and 1.0000   eff. task min and max   -1   20   number of fp tasks:      5
2025-11-13 01:29:23,984 - INFO : -------------------------iter.000000 task 07--------------------------
2025-11-13 01:33:56,374 - INFO : -------------------------iter.000000 task 08--------------------------
2025-11-13 01:33:56,634 - INFO : failed frame:      1 in      5   20.00 %
2025-11-13 01:33:56,660 - INFO : failed tasks:      1 in      5   20.00 %
```

<span style="font-size:16px; color:lightblue;">这里的`thresholds: 1.0000 and 1.0000`是什么意思呢？

0. 第一个1是`fp_accurate_soft_threshold`，第二个1是`fp_accurate_threshold`
1. `fp_accurate_threshold`: 如果准确比率`accurate_ratio`大于fp_accurate_threshold，则不执行02.fp 中fp计算，即fp_task_max = 0
2. `fp_accurate_soft_threshold`: 如果准确率`accurate_ratio`在`fp_accurate_soft_threshold`与`fp_accurate_threshold`之间，那么fp_task_max线性衰减为零。如果准确率`accurate_ratio`小于`fp_accurate_soft_threshold`，满额提交左右的FP任务。


### <span style="color:yellow"> 有时候需要用lammps跑一些特殊的功能，比如蒙卡交换金属原子位置模拟高熵合金的无序性。

这个教程教会我们如何使用model_devi_jobs的template函数提供lammps的模板：https://bohrium-doc.dp.tech/docs/software/DP-GEN_lmp_template/

### <span style="color:yellow"> dpgen通过metadynamics增强采样

1. https://bohrium.dp.tech/notebook/ce0e23b7a5714043824b6e192964b0a5
2. https://zhuanlan.zhihu.com/p/350009924
3. https://bohrium-doc.dp.tech/docs/software/DP-GEN_lmp_template/

需要在param.json中添加开关：
```json
{
"model_devi_plumed": true,
"model_devi_jobs": [
  {"sys_idx":[0], "trj_freq": 10, "_idx": "00",
    "rev_mat":  {"lmp": { "V_NSTEPS": [20000], "V_TEMP": [300], "V_PRES": [700] }, "plm": { "V_TEMP": [300], "V_STRIDE": [10] } },
    "template": {"lmp": "lmp/input.lammps",  "plm": "lmp/input.plumed"}
  }]
}

```