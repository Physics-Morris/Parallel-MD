import h5py
import numpy as np
from pylab import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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

colors = cm.get_cmap('Set1', 8)

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

# ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
# ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
# ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
# ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

def gaus(x, a, x0, sigma):
    return a * np.exp(-(x-x0)**2/(2*sigma**2))

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

vx = particle_velocity[:, 0]
vy = particle_velocity[:, 1]
vz = particle_velocity[:, 2]

hist1 = ax.hist(vx, bins=500, histtype='step', density=True, color=colors(0), label='$v_x$')
hist2 = ax.hist(vy, bins=500, histtype='step', density=True, color=colors(1), label='$v_y$')
hist3 = ax.hist(vz, bins=500, histtype='step', density=True, color=colors(2), label='$v_z$')

# fitting gaussian
# x1, y1 = np.array(hist1[1][:-1]), np.array(hist1[0])
# x2, y2 = np.array(hist2[1][:-1]), np.array(hist2[0])
# x3, y3 = np.array(hist3[1][:-1]), np.array(hist3[0])

# popt1, pcov1 = curve_fit(gaus, x1, y1)
# popt2, pcov2 = curve_fit(gaus, x2, y2)
# popt3, pcov3 = curve_fit(gaus, x3, y3)

# ax.plot(x1, gaus(x1, *popt1), color=colors(3), ls='-', label='$fitting: \sigma$='+str("{:.1f}".format(popt1[2])))
# ax.plot(x2, gaus(x2, *popt2), color=colors(4), ls='-', label='$fitting: \sigma$='+str("{:.1f}".format(popt2[2])))
# ax.plot(x3, gaus(x3, *popt3), color=colors(5), ls='-', label='$fitting: \sigma$='+str("{:.1f}".format(popt3[2])))

ax.legend()
ax.set_xlabel('$v_x, v_y, v_z$')
ax.set_ylabel('$f(v)$')
plt.tight_layout()

# plt.savefig('figures/maxwell_distribution_vx_vy_vz_light.png', dpi=1000)
plt.savefig('figures/maxwell_distribution_vx_vy_vz_dark.png', dpi=1000, transparent=True)