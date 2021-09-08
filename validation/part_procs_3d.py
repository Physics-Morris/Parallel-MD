import h5py
import numpy as np
from pylab import cm
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['figure.figsize'] = [5.6, 4]
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.markersize'] = 6
plt.rcParams['legend.fontsize'] = 13
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.linewidth'] = 1

plt.style.use('dark_background')

colors = cm.get_cmap('tab20b', 27)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='3d')

ax.xaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', top='on')
ax.xaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', top='on')
ax.yaxis.set_tick_params(which='major', size=5, width=1,
                         direction='in', right='on')
ax.yaxis.set_tick_params(which='minor', size=3, width=1,
                         direction='in', right='on')

ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

hdf_start = 0
hdf_end = 0
file_dir = '../data/'
for i in range(hdf_start, hdf_end+1):
    # open file
    hf = h5py.File(file_dir + f'{i:04}.h5', 'r')
    # avaliable keys
    print("Avaliable Keys:", hf.keys())
    # get data
    particle_id = np.array(hf.get('particle_id')[:])
    particle_charge = np.array(hf.get('particle_charge')[:])
    particle_mass = np.array(hf.get('particle_mass')[:])
    particle_position = np.array(hf.get('particle_position')[:])
    particle_velocity = np.array(hf.get('particle_velocity')[:])
    particle_index = np.array(hf.get('global_cell_index')[:])
    particle_procs = np.array(hf.get('processor_rank')[:])

for i in range(particle_id.shape[0]):
    ax.plot(particle_position[i, 0], particle_position[i, 1], particle_position[i, 2], '.', \
            color=colors(particle_procs[i, 0]), markersize=.5, alpha=0.9)

# ax.legend()
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_ylabel('$z$')
plt.tight_layout()
# plt.savefig('figures/uniform_procs_layout_3d_light.png', dpi=1000)
plt.savefig('figures/uniform_procs_layout_3d_dark.png', dpi=1000, transparent=True)