import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle
import scipy
from scipy.stats import gaussian_kde

with open("../../../data/supplementary_experiments/hash_uptake_imaging/hash_intensity_df.pkl", 'rb') as f:
    df = pickle.load(f)

## Supplementary Figure 16a
fig, axis = plt.subplots()
fig.set_size_inches(9,6)
ax = sns.histplot(x='mean_intensity', data = df, bins=50, stat = "density")
ax.set_title('')
ax.set_xlabel('Mean Intensity', fontsize=18)
ax.set_ylabel('Percentage of cells', fontsize=18)

xmin, xmax = plt.xlim()
mu, std = scipy.stats.norm.fit(df['mean_intensity'])
x = np.linspace(xmin, xmax, 100)
p = scipy.stats.norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
plt.show()
fig.savefig("mean_intensity_hist.pdf", bbox_inches='tight')


## Supplementary Figure 16b
x = np.array(df['nuclear_area'].tolist())
y = np.array(df['mean_intensity'].tolist())
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

m, b = np.polyfit(x, y, 1)

fig, ax = plt.subplots()
fig.set_size_inches(9,6)
ax.scatter(x, y, c=z, s=10)
plt.show()

