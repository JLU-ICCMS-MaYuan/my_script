import numpy as np
import matplotlib.pyplot as plt

def read_COHP(fn):
    raw = open(fn).readlines()
    raw = [l for l in raw if 'No' not in l][3:]
    raw = [[eval(i) for i in l.split()] for l in raw]
    return np.array(raw)


data_cohp = read_COHP('COHPCAR.lobster')
labels_cohp = [l[:-1] for l in open('labels').readlines()]
icohp_ef = [eval(l.split()[-1])
            for l in open('ICOHPLIST.lobster').readlines()[1:]]

for i in range(int((len(data_cohp[0])-3)/2)):
    fig, ax1 = plt.subplots(figsize=[2.4, 4.8])
    ax1.plot(-data_cohp[:, i*2+3], data_cohp[:, 0],
             color='k', label=labels_cohp[i])
    ax1.fill_betweenx(data_cohp[:, 0], -data_cohp[:, i*2+3], 0,
                      where=-data_cohp[:, i*2+3] >= 0, facecolor='green', alpha=0.2)
    ax1.fill_betweenx(data_cohp[:, 0], -data_cohp[:, i*2+3], 0,
                      where=-data_cohp[:, i*2+3] <= 0, facecolor='red', alpha=0.2)
    ax1.set_ylim([-10, 6])
    ax1.set_xlim([-0.75, 0.75])
    ax1.set_xlabel('-COHP (eV)', color='k', fontsize='large')
    ax1.set_ylabel('$E-E_\mathrm{F}$ (eV)', fontsize='large')
    ax1.tick_params(axis='x', colors="k")
    # ICOHP
    ax2 = ax1.twiny()
    ax2.plot(-data_cohp[:, i*2+4], data_cohp[:, 0], color='grey')
    ax2.set_ylim([-10, 6])
    ax2.set_xlim([-0.01, 1.5])
    ax2.set_xlabel('-ICOHP (eV)', color='grey', fontsize='large')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.tick_params(axis='x', colors="grey")
    # markers
    ax1.axvline(0, color='k', linestyle=':', alpha=0.5)
    ax1.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax2.annotate(labels_cohp[i], xy=(1.45, 5.5), ha='right', va='top', bbox=dict(
        boxstyle='round', fc='w', alpha=0.5))
    ax2.annotate(f'{-icohp_ef[i]:.3f}', xy=(1.45, -0.05),
                 ha='right', va='top', color='grey')
    fig.savefig(f'cohp-{i+1}.png', dpi=500,
                bbox_inches="tight", transparent=True)
    plt.close()