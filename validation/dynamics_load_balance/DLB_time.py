import numpy as np
from pylab import cm
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['font.size'] = 16
plt.rcParams['figure.figsize'] = [7.6, 4]
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

ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.10))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(.2))

extent = [0.9, 0.5, 0.3, 0.1]
speedup1 = [2.94, 2.79, 2.86, 2.43]
time1 = [2.5, 2.3, 2.1, 2.2]
speedup2 = [2.1, 2.06, 2.04, 1.92]
time2 = [2.1, 1.9, 1.9, 1.8]
speedup3 = [2.74, 2.67, 2.68, 2.38]
time3 = [2.3, 2.0, 2.0, 1.9]
speedup4 = [2.92, 2.70, 2.72, 2.30]
time4 = [2.6, 2.4, 2.5, 2.6]
speedup5 = [3.15, 3.34, 3.21, 2.88]
time5 = [2.3, 2.0, 1.9, 1.8]
speedup6 = [6.86, 6.15, 6.08, 4.74]
time6 = [4.1, 3.6, 3.6, 3.9]
speedup7 = [6.02, 5.96, 5.66, 4.87]
time7 = [3.8, 3.6, 3.4, 3.7]

ax.plot(time1, speedup1, '-o', markersize=8, color=colors(0), label='uniform_slab')
ax.plot(time2, speedup2, '-*', markersize=8, color=colors(1), label='uniform_slab $+$ rotate')
ax.plot(time3, speedup3, '-v', markersize=8, color=colors(2), label='gaussian_x')
ax.plot(time4, speedup4, '-s', markersize=8, color=colors(3), label='gaussian_x $+$ rotate')
ax.plot(time5, speedup5, '-p', markersize=8, color=colors(5), label='gaussian_xy')
ax.plot(time6, speedup6, '-h', markersize=8, color=colors(6), label='gaussian_xy $+$ rotate')
ax.plot(time7, speedup7, '-X', markersize=8, color=colors(7), label='uniform_sphere')

ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel('Time $[s]$')
ax.set_ylabel(r'maximum speedup ($\times$)')
ax.set_title('$10^7$ particles $+$ $8$ processors')

plt.tight_layout()
plt.savefig('time_speedup_dark.png', dpi=1000, transparent=True)
# plt.savefig('time_speedup_light.png', dpi=1000)
plt.show()