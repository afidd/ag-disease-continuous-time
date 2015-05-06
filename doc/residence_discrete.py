'''
This script reads event data from a discrete time simulation
and generates histograms of the number of days for each
transition.
'''
import csv
import logging
import h5py
import numpy as np
import matplotlib.pyplot as plt
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)



def first_dataset(filename, functor):
    f=h5py.File(filename, "r")
    sizes=list()
    infections=set([0, 5, 6])
    for trajgroup in f["/trajectory"]:
        if trajgroup.startswith("dset"):
            functor(trajgroup)
            f.close()
            return
    f.close()


def foreach_dataset(filename, functor):
    logger.debug("Opening {0}".format(filename))
    f=h5py.File(filename, "r")
    sizes=list()
    infections=set([0, 5, 6])
    for trajgroup in f["/trajectory"]:
        if trajgroup.startswith("dset"):
            functor(f["/trajectory/{0}".format(trajgroup)])
    f.close()


class BaseCounts(object):
    def __init__(self):
        self.farm_cnt=-1
        self.run_cnt=0
        self.day_cnt=-1

    def __call__(self, trajgroup):
        farm_cnt=max(trajgroup["Who"])+1
        self.farm_cnt=max(self.farm_cnt, farm_cnt)
        self.day_cnt=max(self.day_cnt, trajgroup["When"][-1])
        self.run_cnt+=1


def binned(observed, censored, measured, end_day):
    delta_t=measured[:,1]-measured[:,0]
    counts=np.histogram(delta_t[delta_t>=0], bins=range(0, len(observed)+1, 1))[0]
    observed+=counts
    has_start=set(np.where(measured[:,0]>=0)[0])
    no_end=set(np.where(measured[:,1]<0)[0])
    cens_farms=list(has_start & no_end)
    cens_counts=np.histogram(end_day-measured[cens_farms,0],
        bins=range(0, len(censored)+1, 1))[0]
    censored+=cens_counts


class Tracking(object):
    def __init__(self, farm_cnt, run_cnt, day_cnt):
        self.farm_cnt=farm_cnt
        self.run_cnt=run_cnt
        # This is a count of times it took this many days.
        self.infect=np.zeros((day_cnt+1,), np.int)
        # The c version is censors.
        self.infectc=np.zeros((day_cnt+1,), np.int)
        self.latent=np.zeros((day_cnt+1,), np.int)
        self.latentc=np.zeros((day_cnt+1,), np.int)
        self.clinical=np.zeros((day_cnt+1,), np.int)
        self.clinicalc=np.zeros((day_cnt+1,), np.int)
        self.run_idx=0

    def __call__(self, trajgroup):
        # each state array is (enabling time, firing time)
        # for each farm. Initialize with -1 when there is no enabling
        # time. Make firing time -2 so that -2 - (-1)= -1, so clearly
        # not a time difference.
        infect=np.ones((self.farm_cnt, 2), np.int)
        infect[:,0]*=0
        infect[:,1]*=-1
        infect[20,0]=-1
        latent=np.ones((self.farm_cnt, 2), np.int)
        latent[:,0]*=-1
        latent[:,1]*=-2
        latent[20,0]=0
        clinical=np.ones((self.farm_cnt, 2), np.int)
        clinical[:,0]=-1
        clinical[:,1]=-2
        events=trajgroup["Event"]
        who=trajgroup["Whom"]
        when=trajgroup["When"]
        transitions_type={(0,1) : 0, (1,3) : 1, (3,4) : 3,
            (4,0) : 4, (0,3) : 5, (0,4) : 6, (1,4) : 8}
        # (0,1), (0,3), (1,3), (3,4)
        today=-1
        for idx in range(len(events)):
            e=events[idx]
            today=when[idx]
            if e==0:
                if (who[idx]==20):
                    logger.debug("20 had an infection")
                if (who[idx]<0):
                    logger.error("writing to a negative who")
                if (who[idx]>=self.farm_cnt):
                    logger.error("farm index too large")
                infect[who[idx]][1]=int(when[idx])
                latent[who[idx]][0]=int(when[idx])
            elif e==5:
                infect[who[idx]][1]=int(when[idx])
                latent[who[idx]][0]=int(when[idx])
                latent[who[idx]][1]=int(when[idx])
                clinical[who[idx]][0]=int(when[idx])
            elif (e==1 or e==7):
                latent[who[idx]][1]=int(when[idx])
                clinical[who[idx]][0]=int(when[idx])
            elif e==3:
                clinical[who[idx]][0]=int(when[idx])
            else:
                logger.error("Unexpected event {0}".format(e))
                assert(False)

        binned(self.infect, self.infectc, infect, today)
        binned(self.latent, self.latentc, latent, today)
        binned(self.clinical, self.clinical, clinical, today)

        self.run_idx+=1


def write_csv(name, observed, censored):
    logger.debug("{0} observed".format(observed))
    logger.debug("{0} censored".format(censored))
    with open("{0}.csv".format(name), "w") as csvfile:
        writer=csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["trial", "value", "censored"])
        idx=1
        for oidx in range(len(observed)):
            for i in range(observed[oidx]):
                writer.writerow([idx, oidx, 1])
                idx+=1
        for oidx in range(len(censored)):
            for i in range(censored[oidx]):
                writer.writerow([idx, oidx, 0])
                idx+=1



if __name__ == "__main__":
    parser=DefaultArgumentParser(description="Finds residence time in states.")
    parser.add_argument("--input", dest="infile", action="store",
        default="naadsm.h5", help="Input HDF5 file with ensemble of events")

    args=parser.parse_args()

    counts=BaseCounts()
    foreach_dataset(args.infile, counts)
    logger.info("Number of farms {0}.".format(counts.farm_cnt))
    logger.info("Number of runs {0}.".format(counts.run_cnt))
    logger.info("Largest number of days {0}.".format(counts.day_cnt))

    tracking=Tracking(counts.farm_cnt, counts.run_cnt, counts.day_cnt)
    foreach_dataset(args.infile, tracking)

    write_csv("susceptible", tracking.infect, tracking.infectc)
    write_csv("latent", tracking.latent, tracking.latentc)
    write_csv("clinical", tracking.clinical, tracking.clinicalc)
