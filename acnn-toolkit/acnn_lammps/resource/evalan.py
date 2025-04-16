#!/usr/bin/env python3
import subprocess
from matplotlib import pyplot
import numpy as np
from sys import argv

#pyplot.rc('font', family='serif', serif='Times New Roman')
pyplot.rcParams['mathtext.fontset'] = 'stix'
pyplot.rcParams["font.family"] = "DejaVu Serif"


def rmse(diff: np.ndarray) -> float:
    """Calculate root mean square error.

    Parameters
    ----------
    diff : np.ndarray
        difference

    Returns
    -------
    float
        root mean square error
    """
    return np.sqrt(np.average(diff * diff))


def get_data(log_file_: str):
    file = open(log_file_)

    lines = file.readlines()

    e_line = None
    e_pred_line = None
    natoms_line = None
    f_line = None
    f_pred_line = None
    v_line = None
    v_pred_line = None
    rmseE = None
    rmseF = None
    rmseV = None

    v = None
    v_pred = None

    for i in range(len(lines)):
        if lines[i] == 'e\n':
            e_line = lines[i + 1]

        if lines[i] == 'e_pred\n':
            e_pred_line = lines[i + 1]

        if lines[i] == 'natoms\n':
            natoms_line = lines[i + 1]

        if lines[i] == 'f:\n':
            f_line = lines[i + 1]

        if lines[i] == 'f_pred:\n':
            f_pred_line = lines[i + 1]

        if 'v:' in lines[i]:
            v_line = lines[i + 1]

        if 'v_pred:' in lines[i]:
            v_pred_line = lines[i + 1]

        if "rmseE" in lines[i]:
            rmseE = float(lines[i].split()[-1])

        if "rmseF" in lines[i]:
            rmseF = float(lines[i].split()[-1])

        if "rmseV" in lines[i]:
            rmseV = float(lines[i].split()[-1])

    e = np.array(e_line.split('[')[1].split(']')[0].split(','), float)
    e_pred = np.array(e_pred_line.split('[')[1].split(']')[0].split(','), float)
    natoms = np.array(natoms_line.split('[')[1].split(']')[0].split(','), float)
    f = np.array(f_line.split('[')[1].split(']')[0].split(','), float)
    f_pred = np.array(f_pred_line.split('[')[1].split(']')[0].split(','), float)
    if v_line:
        v = np.array(v_line.split('[')[1].split(']')[0].split(','), float)
        v_pred = np.array(v_pred_line.split('[')[1].split(']')[0].split(','), float)
    print(e.shape)
    print(e_pred.shape)
    print(natoms.shape)
    print(f.shape)
    print(f_pred.shape)
    if v is not None:
        print(v.shape)
        print(v_pred.shape)
    print(f"rmseE :{rmseE: 16.4e}   eV")
    print(f"rmseF :{rmseF: 16.4e}   eV/A")
    if rmseV is not None:
        print(f"rmseV :{rmseV: 16.4e}   eV")
    return e, e_pred, f, f_pred, v, v_pred, natoms


def plot_fig(e, e_pred, f, f_pred, v, v_pred, natoms):
    fig = pyplot.figure(figsize=(16, 5), dpi=500)

    # Plot Energy
    pyplot.subplot(1, 3, 1)
    pyplot.title('Energy')
    minmin = min(min(e / natoms), min(e_pred / natoms))
    maxmax = max(max(e / natoms), max(e_pred / natoms))
    pyplot.plot([minmin, maxmax], [minmin, maxmax], 'black', linewidth=1)
    pyplot.scatter(e / natoms, e_pred / natoms, 8)
    pyplot.text(minmin + (maxmax - minmin) * 0.5, minmin + (maxmax - minmin) * 0.1,
                f'RMSE\n{rmse((e - e_pred) / natoms):.4e} eV/atom', multialignment='right', fontsize=12)
    pyplot.xlabel("energy by DTF    eV/atom")
    pyplot.ylabel("energy by NN     eV/atom")

    # Plot Force
    pyplot.subplot(1, 3, 2)
    pyplot.title('Force')
    minmin = min(min(f), min(f_pred))
    maxmax = max(max(f), max(f_pred))
    pyplot.plot([minmin, maxmax], [minmin, maxmax], 'black', linewidth=1)
    pyplot.scatter(f, f_pred, 8, 'm')
    pyplot.text(minmin + (maxmax - minmin) * 0.5, minmin + (maxmax - minmin) * 0.1,
                f'RMSE\n{rmse(f - f_pred):.4e} eV/ $\AA$', multialignment='right', fontsize=12)
    pyplot.xlabel(r"force by DTF    eV/$\AA$")
    pyplot.ylabel(r"force by NN     eV/$\AA$")

    # Plot Stress
    if v is not None and v_pred is not None:
        pyplot.subplot(1, 3, 3)
        pyplot.title('Stress')
        minmin = min(min(v / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1)), min(v_pred / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1)))
        maxmax = max(max(v / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1)), max(v_pred / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1)))
        pyplot.plot([minmin, maxmax], [minmin, maxmax], 'black', linewidth=1)
        pyplot.scatter(v / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1), v_pred / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1), 8, 'g')
        pyplot.text(minmin + (maxmax - minmin) * 0.5, minmin + (maxmax - minmin) * 0.1,
                    f'RMSE\n{rmse((v - v_pred) / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1)):.4e} eV', multialignment='right', fontsize=12)
        pyplot.xlabel("stress by DTF    eV")
        pyplot.ylabel("stress by NN     eV")

    pyplot.savefig('plot')


def get_xsf(log_file_):
    grep_command = f"grep xsf {log_file_}"
    xsf_file = subprocess.check_output(grep_command, shell=True, text=True).split()
    return xsf_file


def get_f_ends(i_, natoms_):
    if i_ == 0:
        return 0, int(natoms_[0] * 3)

    cum_natoms = np.cumsum(natoms_)
    return int(cum_natoms[i_ - 1] * 3), int(cum_natoms[i_] * 3)


def get_v_ends(i_):
    return int(i_ * 9), int((i_ + 1) * 9)


if __name__ == '__main__':
    log_file = argv[1]
    print(f"analysis {log_file} ...")
    e, e_pred, f, f_pred, v, v_pred, natoms = get_data(log_file)

    rmseE = rmse((e - e_pred) / natoms)
    rmseF = rmse(f - f_pred)
    print(f"rmseE:{rmseE: 16.4e}   eV")
    print(f"rmseF:{rmseF: 16.4e}   eV/A")
    if v is not None:
        rmseV = rmse((v - v_pred) / np.broadcast_to(natoms, [9, natoms.size]).T.reshape(-1))
        print(f"rmseV:{rmseV: 16.4e}   eV")

    plot_fig(e, e_pred, f, f_pred, v, v_pred, natoms)

    xsf_files = get_xsf(log_file)
    print(" " * 72, "e/atom      f/atom      v/atom")
    for i in range(len(xsf_files)):
        be, ee = get_f_ends(i, natoms)
        bv, fv = get_v_ends(i)
        print(
            f"{xsf_files[i]:64}  {int(natoms[i]):d}  {rmse((e[i] - e_pred[i]) / natoms[i]):.4e}  {rmse(f[be:ee] - f_pred[be:ee]):.4e}", end="")
        if v is not None:
            print(f"  {rmse((v[bv:fv] - v_pred[bv:fv]) / natoms[i]):.4e}")
        else:
            print()

