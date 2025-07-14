import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('E://Full//MSS t (0-129, 4w).csv')

plt.figure(figsize=(5, 2.5))
plt.plot(data['Pos. [°2θ]'], data['Intensity [Counts]'],lw=.25)
plt.yscale('log')
plt.ylim(300, max(data['Intensity [Counts]'][200:]) * 1.1)  # Set y-axis limit to 10% above max intensity
plt.xlim(6, max(data['Pos. [°2θ]']) + 1)  # Set x-axis limit to 1 degree above max position
# plt.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
plt.title('XRD Scan')
plt.xlabel('2θ Position (°)')
plt.ylabel('Intensity (Counts)')
plt.show()