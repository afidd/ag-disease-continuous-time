"""
Rain simulation

Simulates rain drops on a surface by animating the scale and opacity
of 50 scatter points.

Author: Nicolas P. Rougier
"""
import logging
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from  matplotlib.animation import FuncAnimation


def datasets(open_file):
    trajectories=open_file['/trajectory']
    names=list()
    for trajectory_name in trajectories:
        if trajectory_name.startswith("dset"):
            names.append(trajectory_name)
    return names


def transitions(open_file):
    locations=open_file["/trajectory/locations"]
    trajectory_name=datasets(open_file)[0]
    dset=open_file["/trajectory/{0}".format(trajectory_name)]
    events=list()
    times=list()
    events.append(0)
    times.append(0.0)
    for kind, whom, who, when in dset:
        if kind==0:
            events.append(whom)
            times.append(when)
    return (locations, np.array(events), np.array(times))


logger=logging.getLogger(__file__)
logging.basicConfig(level=logging.DEBUG)

data_file="src/sirexp.h5"
f=h5py.File(data_file)
locations, farm_order, farm_times=transitions(f)

color_choice={'susceptible' : 'black', 'infected' : 'orangered'}
color_code=dict()
for disease_name, cname in color_choice.items():
    r, g, b=mcolors.hex2color(mcolors.cnames[cname])
    color_code[disease_name]=np.array((r, g, b, 1.0))

# Create new Figure and an Axes which fills it.
fig = plt.figure(figsize=(7,7))
ax = fig.add_axes([0, 0, 1, 1], frameon=False)
ax.set_xlim(0,1), ax.set_xticks([])
ax.set_ylim(0,1), ax.set_yticks([])

# Create rain data
n_drops = locations.shape[0]
farm_cnt = n_drops
end_time=farm_times[-1]*1.01
frame_cnt=1000
frame_interval=10

# Create disease data
#farm_order=np.arange(0, farm_cnt)
#np.random.shuffle(farm_order)
#farm_times=np.random.uniform(0, end_time, farm_cnt)
#farm_times.sort()
current_farm=0


farms = np.zeros(n_drops, dtype=[('position', float, 2),
                                      ('size',     float, 1),
                                      ('color',    float, 4)])

marker_size=100

# Initialize the raindrops in random positions and with
# random growth rates.
farms['position'] = locations #np.random.uniform(0, 1, (n_drops, 2))
farms['size'].fill(marker_size)
farms['color'][:]=color_code['susceptible']

rain_growth_rate=50

rain_drops = np.zeros(n_drops, dtype=[('position', float, 2),
                                      ('size',     float, 1),
                                      ('color',    float, 4)])

rain_drops['position'] = farms['position']
rain_drops['size'].fill(marker_size)
rain_drops['color'][:]=np.array((0.0, 0.0, 0.0, 0.0))



# Construct the scatter which we will update during animation
# as the raindrops develop.
farms_scat = ax.scatter(farms['position'][:,0], farms['position'][:,1],
                  s=farms['size'], lw=0.5, facecolors=farms['color'],
                  edgecolors='none')
rain_scat = ax.scatter(rain_drops['position'][:,0], rain_drops['position'][:,1],
                  s=rain_drops['size'], lw=0.5, facecolors='none',
                  edgecolors=rain_drops['color'])


def update(frame_number):
    # Get an index which we can use to re-spawn the oldest raindrop.
    current_time = frame_number*end_time/frame_cnt

    # Make all colors more transparent as time progresses.
    rain_drops['color'][:, 3] -= 1.0/len(rain_drops)
    rain_drops['color'][:,3] = np.clip(rain_drops['color'][:,3], 0, 1)

    # Make all circles bigger.
    rain_drops['size'] += rain_growth_rate

    global current_farm
    global farm_order
    global farm_times
    global farm_cnt
    while current_farm<len(farm_order) and \
            farm_times[current_farm]<current_time:
        farm_idx=farm_order[current_farm]
        logger.debug("changing farm {0}".format(farm_idx))
        farms['color'][farm_idx] = color_code['infected']
        rain_drops['size'][farm_idx]=marker_size
        rain_drops['color'][farm_idx, 3]=1.0
        current_farm+=1


    # Update the scatter collection, with the new colors, sizes and positions.
    farms_scat.set_facecolors(farms['color'])
    rain_scat.set_edgecolors(rain_drops['color'])
    rain_scat.set_sizes(rain_drops['size'])


# Construct the animation, using the update function as the animation
# director.
animation = FuncAnimation(fig, update, frames=frame_cnt,
    interval=frame_interval)
animation.save("points.mp4")
#plt.show()
