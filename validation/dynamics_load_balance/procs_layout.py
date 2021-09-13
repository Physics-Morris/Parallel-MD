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
mpl.rcParams['figure.dpi'] = 150

# plt.style.use('dark_background')

colors = cm.get_cmap('tab20', 10)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

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

hdf_start = 1
hdf_end = 1
file_dir = 'data/'
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
    particle_index = np.array(hf.get('particle_global_cell_index')[:])
    particle_procs = np.array(hf.get('particle_processor_rank')[:])

x_slice = np.array(hf.get('simulation_info/DLB_slice/x_slice'))
y_slice = np.array(hf.get('simulation_info/DLB_slice/y_slice'))
z_slice = np.array(hf.get('simulation_info/DLB_slice/z_slice'))
cell_num = np.array(hf.get('simulation_info/DLB_slice/cell_number'))

for i in range(particle_id.shape[0]):
    ax.plot(particle_position[i, 0], particle_position[i, 1], '.', \
            color=colors(particle_procs[i, 0]), markersize=0.8, alpha=0.8)

# plot processor layout
cell_wx = 100 / cell_num[0]
cell_wy = 100 / cell_num[1]
cell_wz = 100 / cell_num[2]

# draw y slice
ax.hlines((y_slice[0][0]-1)*cell_wy, 20, 70, ls='--', color='red', lw=1.5)
ax.hlines((y_slice[1][0]-1)*cell_wy, 20, 70, ls='--', color='red', lw=1.5)
ax.hlines((y_slice[2][0]-1)*cell_wy, 20, 70, ls='--', color='red', lw=1.5)
ax.hlines((y_slice[3][0]-1)*cell_wy, 20, 70, ls='--', color='red', lw=1.5)
for i in range(1, cell_num[1][0]):
    ax.hlines(i*cell_wy, 20, 70, ls='-', color='gray', alpha=0.3, lw=0.5)

# draw x slice
ax.vlines((x_slice[0][0][0]-1)*cell_wx, 20, (y_slice[0][0]-1)*cell_wy, ls='--', color='red', lw=1.5)
ax.vlines((x_slice[0][1][0]-1)*cell_wx, (y_slice[0][0]-1)*cell_wy, (y_slice[1][0]-1)*cell_wy, ls='--', color='red', lw=1.5)
ax.vlines((x_slice[0][2][0]-1)*cell_wx, (y_slice[1][0]-1)*cell_wy, (y_slice[2][0]-1)*cell_wy, ls='--', color='red', lw=1.5)
ax.vlines((x_slice[0][3][0]-1)*cell_wx, (y_slice[2][0]-1)*cell_wy, (y_slice[3][0]-1)*cell_wy, ls='--', color='red', lw=1.5)
ax.vlines((x_slice[0][4][0]-1)*cell_wx, (y_slice[3][0]-1)*cell_wy, 70, ls='--', color='red', lw=1.5)
for i in range(1, cell_num[0][0]):
    ax.vlines(i*cell_wx, 20, 70, ls='-', color='gray', alpha=0.3, lw=0.5)

ax.set_xlim(20, 70)
ax.set_ylim(20, 70)
ax.set_title(r'$[2 \times 5 \times 1]$, $dlb\_extent=0.1$')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_aspect('equal')
plt.tight_layout()
# plt.savefig('procs_layout2_01_light.png', dpi=1000)
# plt.savefig('procs_layout2_01_dark.png', dpi=1000, transparent=True)
plt.show()