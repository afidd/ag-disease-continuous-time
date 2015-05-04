'''
This plots total outbreak size for an ensemble.
Working from the pretty plot at
http://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d.html
'''
import csv
import h5py
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

def run_sizes(filename):
    f=h5py.File(filename)
    sizes=list()
    for trajgroup in f["/trajectory"]:
        if trajgroup.startswith("dset"):
            cnt=len(f["/trajectory/{0}/Event".format(trajgroup)])
            sizes.append(cnt//3)
    return sizes

totals=run_sizes("../run.h5")
with open('sizesc.csv', 'w') as csvfile:
    writer=csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
    writer.writerow(["trial", "outbreaksize"])
    for i in range(len(totals)):
        writer.writerow([i+1, totals[i]])

X_plot=np.linspace(-5, 50, 1000)[:, np.newaxis]
fig, ax=plt.subplots(1, 1)

# Gaussian KDE
kde=KernelDensity(kernel="gaussian", bandwidth=3).fit(totals)
log_dens = kde.score_samples(X_plot)
ax[0,0].fill(X_plot[:, 0], np.exp(log_dens), fc='#AAAFF')

