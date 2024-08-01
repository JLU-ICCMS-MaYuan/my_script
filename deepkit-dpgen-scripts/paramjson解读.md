```shell
{
	"model_devi_engine":"calypso",
	"calypso_input_path":"/public/home/mayuan/phycollege_workplace/5.calypso/22.Sc-Ca-H/0.MLP/calypso_input",
	"model_devi_max_iter": 15,
	"vsc":true,
	"ratio_failed":0.2,
	"type_map": ["Sc", "Ca","H"],
	"mass_map": [45, 40, 1],
	"init_data_prefix":"22.Sc-Ca-H/0.MLP/1.dp-rawdata", # init_data_prefix 是dp格式数据集的路径的前缀
													    # init_data_sys 是dp格式数据集的路径下的具体各个数据的路径
	"init_data_sys": [ 
		"Sc1Ca1H1",
		"Sc1Ca1H10",
		"Sc1Ca1H11",
		"Sc1Ca1H12",
		"Sc1Ca1H13",
		"Sc1Ca1H14",
		"Sc5Ca5H"
    ],
    "_comment": " that's all ",
	"_comment3": "sys_configs_prefix 指定跑分子动力学的初始结构的路径的前缀。sys_configs在上述路径下的具体结构的路径",
	"sys_configs_prefix": "",
    "sys_configs": "", 

    "training_init_model": true,
    "training_iter0_model_path": "/public/home/mayuan/phycollege_workplace/5.calypso/22.Sc-Ca-H/0.MLP/3.adjust/iter.000000/00[0-4]",
    "training_reuse_iter":              1,
    "training_reuse_old_ratio":         0.9,
    "training_reuse_start_lr":          1e-4,
    "training_reuse_stop_batch":        1000000,
    "training_reuse_start_pref_e":      0.2,
    "training_reuse_start_pref_f":      100,

    "numb_models": 4,
    "default_training_param": {
    "model": {
       	"_use_srtab": "/home/wangzy/workplace/dpgen/DP-BJ-PBS/LiLaH-adjust/ZBL/ZBL-Tab",
        "_use_srtab": "../../../ZBL/ZBL-Tab",
       	"_smin_alpha":   0.1,
       	"_sw_rmin":      0.58,
       	"_sw_rmax":      0.68, 
        "descriptor": {
            "type": "se_e2_a",
            "sel": [
                1100,
                1100,
                1500
            ],
            "rcut_smth": 2.0,
            "rcut": 6.0,
            "neuron": [
                25,
                50,
                100
            ],
            "resnet_dt": false,
            "axis_neuron": 12,
            "type_one_side": true,
            "seed": 1801819940,
            "_activation_function": "tanh"
        },
        "fitting_net": {
            "neuron": [
                240,
                240,
                240
            ],
            "resnet_dt": true,
            "seed": 2375417769
        },
        "type_map": [
            "Sc",
            "Ca",
            "H"
        ]
    },
    "learning_rate": {
        "type": "exp",
        "start_lr": 0.001,
        "decay_steps": 5000
    },
    "loss": {
        "start_pref_e": 0.2,
        "limit_pref_e": 2,
        "start_pref_f": 100,
        "limit_pref_f": 1,
        "start_pref_v": 0,
        "limit_pref_v": 0 
    },
    "training": {
        "seed": 3982377700,
        "disp_file": "lcurve.out",
        "disp_freq": 2000,
        "numb_test": 4,
        "save_freq": 2000,
        "save_ckpt": "model.ckpt",
        "disp_training": true,
        "time_training": true,
	"_tensorboard":	true,
	"_tensorboard_log_dir":"log",
	"_tensorboard_freq": 1000,
        "profiling": false,
        "profiling_file": "timeline.json",
        "_set_prefix": "set"
    }
	},
    "model_devi_dt": 0.002,        # 0.002ps=2fs
    "model_devi_skip": 0,
    "model_devi_f_trust_lo": 0.5,  # 挑选结构时力误差最大偏差的下限 eV/A
    "model_devi_f_trust_hi": 0.9,  # 挑选结构时力误差最大偏差的上限 eV/A
    "model_devi_e_trust_lo": 10000000000.0,
    "model_devi_e_trust_hi": 10000000000.0,
    "model_devi_clean_traj": true, # 清除MD的
    "model_devi_jobs": [
    ],

    "fp_style": "vasp",
    "shuffle_poscar": false, # 对POSCAR打乱，无所谓的参数
    "fp_task_max": 10,       # 最多挑fp_task_max个点进行第一性原理计算
    "fp_task_min": 1,        # 最少挑fp_task_min个点进行第一性原理计算
    "fp_pp_path": "/public/home/mayuan/phycollege_workplace/5.calypso/22.Sc-Ca-H/0.MLP/vasp_input",
    "fp_pp_files": [
		"POTCAR.Sc",
        "POTCAR.Ca",
        "POTCAR.H"
    ],
    "fp_incar":"/public/home/mayuan/phycollege_workplace/5.calypso/22.Sc-Ca-H/0.MLP/vasp_input/INCAR"
}
```