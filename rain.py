"""
Rain simulation

Simulates rain drops on a surface by animating the scale and opacity
of 50 scatter points.

Author: Nicolas P. Rougier
"""
import logging
import sys
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

    event=dset["Event"]
    who=dset["Who"]
    whom=dset["Whom"]
    when=dset["When"]
    ret=[np.array(locations)]
    ret.extend([np.array(x) for x in [event, who, whom, when]])
    return ret



_degrees_to_radians=np.pi/180
_radians_km=180*60*1.852/np.pi

def distancekm(latlon1, latlon2):
    """Distance computed on a spherical earth.
    Taken from http://williams.best.vwh.net/avform.htm."""
    ll1=latlon1*_degrees_to_radians
    ll2=latlon2*_degrees_to_radians
    return _radians_km*(2*np.arcsin(np.sqrt(np.power(np.sin((ll1[0]-ll2[0])/2),2)+
        np.cos(ll1[0])*np.cos(ll2[0])*np.power(np.sin((ll1[1]-ll2[1])/2), 2))))


def map_proj():
    import pyproj
    # This is a projection for South Carolina
    sc_proj='+proj=lcc +lat_1=32.5 +lat_2=34.83333333333334 +lat_0=31.83333333333333 +lon_0=-81 +x_0=609600 +y_0=0 +datum=NAD83 +units=ft +no_defs '
    projection=pyproj.Proj(sc_proj)
    return projection

logger=logging.getLogger(__file__)
logging.basicConfig(level=logging.DEBUG)

data_file="run.h5"
f=h5py.File(data_file)
locations, event, who, whom, when=transitions(f)
logger.debug("event {0}".format(event))
logger.debug("who {0}".format(who))
logger.debug("whom {0}".format(whom))
logger.debug("when {0}".format(when))

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
projection=map_proj()
projected=locations[:,:]
for i in range(locations.shape[0]):
    projected[i,:]=np.array(projection(locations[i,1], locations[i,0]))
logger.debug("projected {0}".format(projected))

max_x=np.max(projected[:,0])
min_x=np.min(projected[:,0])
max_y=np.max(projected[:,1])
min_y=np.min(projected[:,1])
locations_scaled=(projected-np.array([min_x, min_y]))*np.array(
    [1.0/(max_x-min_x), 1.0/(max_y-min_y)])
#logger.debug(locations_scaled)
logger.debug("min farm {0} max farm {1}".format(np.min(whom), np.max(whom)))
farm_cnt = n_drops
last_infection=np.where(event==0)[0]
end_time=when[last_infection[-1]]*1.05
frame_cnt=1000
frame_interval=10

# Create disease data
#farm_order=np.arange(0, farm_cnt)
#np.random.shuffle(farm_order)
#farm_times=np.random.uniform(0, end_time, farm_cnt)
#farm_times.sort()
event_idx=0

farms = np.zeros(n_drops, dtype=[('position', float, 2),
                                      ('size',     float, 1),
                                      ('color',    float, 4)])

marker_size=100

# Initialize the raindrops in random positions and with
# random growth rates.
farms['position'] = locations_scaled #np.random.uniform(0, 1, (n_drops, 2))
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

farm_idx=21-1
farms['color'][farm_idx] = color_code['infected']
rain_drops['size'][farm_idx]=marker_size
rain_drops['color'][farm_idx, 3]=1.0

def update(frame_number):
    # Get an index which we can use to re-spawn the oldest raindrop.
    current_time = frame_number*end_time/frame_cnt

    # Make all colors more transparent as time progresses.
    rain_drops['color'][:, 3] -= 1.0/len(rain_drops)
    rain_drops['color'][:,3] = np.clip(rain_drops['color'][:,3], 0, 1)

    # Make all circles bigger.
    rain_drops['size'] += rain_growth_rate
    (einfect, eclinical, erecover, ewane)=(0, 1, 2, 3)

    global event_idx
    global locations_scaled
    global ax
    global event
    global when
    global who
    global whom
    while event_idx<len(event) and \
            when[event_idx]<current_time:
        farm_idx=whom[event_idx]-1
        source_idx=who[event_idx]-1
        now=when[event_idx]
        what=event[event_idx]
        if what==einfect:
            logger.debug("changing farm {0}".format(farm_idx))
            farms['color'][farm_idx] = color_code['infected']
            rain_drops['size'][farm_idx]=marker_size
            rain_drops['color'][farm_idx, 3]=1.0
            x0, y0=locations_scaled[source_idx,:]
            x1, y1=locations_scaled[farm_idx,:]
            dx=x1-x0
            dy=y1-y0
            r=np.sqrt(dx*dx+dy*dy)
            dr=0.04 # Want to reduce length by fixed distance
            g=(r-dr)/r
            if g<0:
                g=0.1
            logger.info("x0 {0} x1 {1}".format(x0, x1))
            ax.arrow(x0, y0, g*dx, g*dy, head_width=0.01, head_length=0.02,
                fc='k', ec='k')
        event_idx+=1


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
