import h5py 
import numpy as np
import os

def compress(params):
    
    if not os.path.exists(params.pos_file):
                    print('\nERROR: file {} not found\n'.format(params.pos_file))
                    exit()

    print('\nCompressing position data into .hdf5 database\n'
            'This may take a while...\n')

    if params.file_format == 'lammps':
        num_steps = (params.num_steps//params.stride) 
        with h5py.File(params.database_file,'w') as fout:
            with open(params.pos_file,'r') as fin: 
                pos_dset = fout.create_dataset('pos',
                        (num_steps,params.num_atoms,3))
                for i in range(num_steps):
                    pos = np.zeros((params.num_atoms,3))
                    for j in range(9):
                        fin.readline()
                    for j in range(params.num_atoms):
                       pos[j,:] = fin.readline().strip().split()[2:]
                    pos_dset[i,:,:] = pos

    print('\nDone compressing {} into .hdf5 format.'
         '\nThe compressed file is \'{}\' (DON\'T CHANGE IT!)\n\n'
         'Set COMPRESS = 0 in the input file and run the code again to continue\n'
         .format(params.pos_file,params.database_file))


