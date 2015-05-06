'''
Produces a csv of outbreak sizes.
'''
import csv
import logging
import h5py
import numpy as np
import matplotlib.pyplot as plt
from default_parser import DefaultArgumentParser
# from sklearn.neighbors import KernelDensity

logger=logging.getLogger(__file__)


def run_sizes(filename):
    f=h5py.File(filename)
    sizes=list()
    infections=set([0, 5, 6])
    for trajgroup in f["/trajectory"]:
        if trajgroup.startswith("dset"):
            cnt=0
            events=f["/trajectory/{0}/Event".format(trajgroup)]
            cnt=0
            for e in events:
                if e in infections:
                    cnt+=1
            sizes.append(cnt)
    return sizes

def write_totals(filename, outfile):
    logger.info("Reading input {0}. Writing to {1}".format(filename, outfile))
    totals=run_sizes(filename)
    logger.info("Number of trajectories {0}, average size {1}".format(
        len(totals), np.average(totals)))
    logger.debug("Sizes are {0}".format(totals))
    with open(outfile, 'w') as csvfile:
        writer=csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["trial", "outbreaksize"])
        for i in range(len(totals)):
            writer.writerow([i+1, totals[i]])

# X_plot=np.linspace(-5, 50, 1000)[:, np.newaxis]
# fig, ax=plt.subplots(1, 1)

# # Gaussian KDE
# kde=KernelDensity(kernel="gaussian", bandwidth=3).fit(totals)
# log_dens = kde.score_samples(X_plot)
# ax[0,0].fill(X_plot[:, 0], np.exp(log_dens), fc='#AAAFF')

if __name__ == '__main__':
    parser=DefaultArgumentParser(description="Produces csv of total outbreak size")
    parser.add_argument("--input", dest="infile", action="store",
        default="run.h5", help="Input HDF5 file with ensemble of events")
    parser.add_argument("--output", dest="outfile", action="store",
        default="sizesc.csv", help="CSV output with sizes")

    args=parser.parse_args()
    write_totals(args.infile, args.outfile)
