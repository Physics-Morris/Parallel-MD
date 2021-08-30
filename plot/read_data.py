import h5py
import numpy as np

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