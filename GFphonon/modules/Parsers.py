import numpy as np
import os
import h5py
import yaml
from yaml import CLoader as Loader

def print_error(txt):
    print('\nERROR: value for input paramater {} seems wrong\n'.format(txt))
    exit()

class parse_input:
    def __init__(self,input_file):
        
        ### defaults
        self.compress = False
        self.debug = False
        self.time_step = 0.5*1e-15
        self.stride = 1
        self.num_splits = 1
        self.num_qpoints = 1
        self.qpoints = np.array([[0.0,0.0,0.0]])
        self.file_format = 'lammps'
        self.out_prefix = 'gf_phonons'
        self.lattice_file = 'lattice.dat'
        self.pos_file = 'pos.dat'
        self.lat_params = [5.431,5.431,5.431]
        self.prim_vecs = np.array(([[1,0,0],[0,1,0],[0,0,1]]))
        self.database_file = 'pos.hdf5'
        self.temperature = 300
        self.plot = False

        self.input_file = input_file
        input_txt = open(input_file,'r').readlines()
        for line in input_txt:
            txt = line.strip()

            # skip blank and comment lines
            if len(line.strip()) == 0: 
                continue
            elif line.strip()[0] == '#':
                continue
            else:
                txt = txt.split()

            # crystal
            if txt[0] == 'LAT_PARAMS':
                try:
                    self.lat_params = np.array(txt[(txt.index('=')+1):
                                            (txt.index('=')+4)]).astype(float)
                except:
                    print_error('LAT_PARAMS')
            elif txt[0] == 'PRIM_VECS':
                try:
                    self.prim_vecs = np.array(txt[(txt.index('=')+1):
                        (txt.index('=')+10)]).astype(float)
                    self.prim_vecs = self.prim_vecs.reshape(3,3)
                except:
                    print_error('PRIM_VECS')
            
            # Q-points
            elif txt[0] == 'NUM_QPOINTS':
                try:
                    self.num_qpoints = int(txt[(txt.index('=')+1)])
                except:
                    print_error('NUM_QPOINTS')
            elif txt[0] == 'QPOINTS':
                try:
                    self.qpoints = np.array(txt[(txt.index('=')+1):
                        (txt.index('=')+int((self.num_qpoints)*3+1))]).astype(float)
                    self.qpoints = self.qpoints.reshape(self.num_qpoints,3)
                except:
                    print_error('QPOINTS')

            # simulation control parameters
            elif txt[0] == 'NUM_STEPS':
                try:
                    self.num_steps = int(txt[txt.index('=')+1])
                except:
                    print_error('NUM_STEPS')
            elif txt[0] == 'TEMPERATURE':
                try:
                    self.temperature = int(txt[txt.index('=')+1])
                except:
                    print_error('TEMPERATURE')
            elif txt[0] == 'STRIDE':
                try:
                    self.stride = int(txt[txt.index('=')+1])
                except:
                    print_error('STRIDE')
            elif txt[0] == 'NUM_ENSEMBLES':
                try:
                    self.num_ensembles = int(txt[txt.index('=')+1])
                except:
                    print_error('NUM_ENSEMBLES')
            elif txt[0] == 'TIME_STEP':
                try:
                    self.time_step = float(txt[txt.index('=')+1])
                except:
                    print_error('TIME_STEP')
                self.time_step = self.time_step*1e-15
            # options
            elif txt[0] == 'COMPRESS':
                try:
                    self.compress = bool(int(txt[txt.index('=')+1]))
                except:
                    print_error('COMPRESS')
            elif txt[0] == 'DEBUG':
                try:
                    self.debug = bool(int(txt[txt.index('=')+1]))
                except:
                    print_error('DEBUG')
            elif txt[0] == 'PLOT':
                try:
                    self.plot = bool(int(txt[txt.index('=')+1]))
                except:
                    print_error('PLOT')                    
            # file names
            elif txt[0] == 'FILE_FORMAT':
                try:
                    self.file_format = str(txt[txt.index('=')+1].strip('\''))
                except:
                    print_error('FILE_FORMAT')
            elif txt[0] == 'POS_FILE':
                try:
                    self.pos_file = str(txt[txt.index('=')+1].strip('\''))
                except:
                    print_error('POS_FILE')
            elif txt[0] == 'LATTICE_FILE':
                try:
                    self.lattice_file = str(txt[txt.index('=')+1].strip('\''))
                except:
                    print_error('LATTICE_FILE')
                if not os.path.exists(self.lattice_file):
                    print('\nERROR: file {} not found\n'.format(self.lattice_file))
                    exit()
            elif txt[0] == 'OUT_PREFIX':
                try:
                    self.out_prefix = str(txt[txt.index('=')+1].strip('\''))
                except:
                    print_error('OUT_PREFIX')
            else: 
                print('\nERROR: option {} not recognized\n'.format(txt[0]))
#                exit()

class parse_lattice:
    def __init__(self,params):
        if not os.path.exists(params.lattice_file):
            print('\nERROR: file {} not found\n'.format(self.lattice_file))
            exit()
        
        atom_ids, unit_cells, basis_pos, masses = np.loadtxt(
                params.lattice_file,unpack=True,skiprows=1)

        with open(params.lattice_file,'r') as fid:
            tmp = fid.readline().strip().split()
            self.num_atoms = int(tmp[0])
            self.num_unitcell = int(tmp[1])
            self.num_basis = int(tmp[2])

        # set proper data types
        self.atom_ids = atom_ids.astype(int)
        self.unit_cells = unit_cells.astype(int)
        self.basis_pos = basis_pos.astype(int)
        self.masses = masses.astype(float)

class parse_eigen_vecs:
    def __init__(self,params):
        if params.with_eigs: # if WITH_EIGS = 1 in input file
            phonopy_data = yaml.load(open(params.eigvecs_file),Loader=Loader)
            
            # data from the phonopy output file
            self.natom = phonopy_data['natom']
            self.num_qpoints = phonopy_data['nqpoint']
            self.qpoints = np.zeros((self.num_qpoints,3))

            # incase the order of basis atoms in phonopy don't match
            if len(params.basis_list) == 0:
                basis_slice = np.arange(self.natom)
            else:
                basis_slice = sorted(params.basis_list)
                for i in range(len(basis_slice)):
                    basis_slice[i] = basis_slice[i]-1
        
            # dispersions from phonopy
            self.freq = np.zeros((self.num_qpoints,self.natom*3))
            self.eig_vecs = np.zeros((self.num_qpoints,self.natom*3,
                self.natom,3)).astype(complex)

            # look up the eigenvectors
            for i in range(self.num_qpoints):
                # loop over q-points
                self.qpoints[i,:] = phonopy_data['phonon'][i]['q-position']

                # loop over bands
                for j in range(self.natom*3):
                    self.freq[i][j] = phonopy_data['phonon'][i]['band'][j]['frequency']

                    # loop over basis atoms
                    for ind in basis_slice:
                            self.eig_vecs[i,j,ind,0] = (
                        phonopy_data['phonon'][i]['band'][j]['eigenvector'][ind][0][0]+
                        1j*phonopy_data['phonon'][i]['band'][j]['eigenvector'][ind][0][1]
                            )
                            self.eig_vecs[i,j,ind,1] = (
                        phonopy_data['phonon'][i]['band'][j]['eigenvector'][ind][1][0]+
                        1j*phonopy_data['phonon'][i]['band'][j]['eigenvector'][ind][1][1]
                            )
                            self.eig_vecs[i,j,ind,2] = (
                        phonopy_data['phonon'][i]['band'][j]['eigenvector'][ind][2][0]+
                        1j*phonopy_data['phonon'][i]['band'][j]['eigenvector'][ind][2][1]
                            )

        else: # intialize empty array 
            self.eig_vecs = []


