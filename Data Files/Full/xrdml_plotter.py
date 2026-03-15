import numpy as np
import matplotlib.pyplot as plt

data45 = np.genfromtxt("Data Files\Full\MSS p (0-129, 3w)_45rotation.txt")
data90 = np.genfromtxt("Data Files\Full\MSS p (0-129, 3w)_90rotation.txt")
data00 = np.genfromtxt("Data Files\Full\MSS p (0-129, 3w).txt")

angles = np.linspace(0.04,128.92,len(data45))

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=True)

ax1.plot(angles, data45, label = "45")
ax1.set_title("45 degree")

ax2.plot(angles, data90, label = "90")
ax2.set_title("90 degree")

ax3.plot(angles,data00, label = "0")
ax3.set_title("00 degree")

fig.supylabel('Intensity (counts)')
fig.supxlabel('Angle (°)')
plt.yscale('log')
plt.tight_layout()
plt.savefig("Comparison_of_Sample_preferred_Orientation_3w.svg")
# plt.show()