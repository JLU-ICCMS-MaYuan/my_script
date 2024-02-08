import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

with open("lcurve.out") as f:
    headers = f.readline().split()[1:]
lcurve = pd.DataFrame(np.loadtxt("lcurve.out"), columns=headers)
legends = ["rmse_e_val", "rmse_e_trn", "rmse_f_val" , "rmse_f_trn" ]
for legend in legends:
    plt.loglog(lcurve["step"], lcurve[legend], label = legend )
    plt.legend()
    plt.xlabel("Training steps")
    plt.ylabel("Loss")
    plt.show()
    plt.savefig('lcurve.png')
