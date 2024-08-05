import sys
import dpdata
import matplotlib.pyplot as plt
import numpy as np 


# 制定训练集的路径
traindata_path = sys.argv[1]
# 制定模型路径
pbfile_path = sys.argv[2]
training_systems = dpdata.LabeledSystem(traindata_path, mt = "deepmd/npy")
predict = training_systems.predict(pbfile_path)

plt.scatter(training_systems['energies'], predict['energies'])

x_range = np.linspace(plt.xlim()[0], plt.xlim()[1])

plt.plot(x_range, x_range, 'r--', linewidth=0.25)
plt.xlabel("Energy of DFT")
plt.ylable("Energy predicted by MPL")
plt.show()
plt.savefig("DFT_vs_MPL.png")
