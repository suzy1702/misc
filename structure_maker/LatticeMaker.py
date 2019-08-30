import numpy as np

def error(string):
    print('\nERROR: {}!\n'.format(string))
    exit()

class lattice:
    def __init__(self,description='Crystal Lattice'):
        self.description = str(description)

    def basis(self,positions,basis_types,masses=False,reduced_coords=True):
        try:
            self.basis_positions = np.array(positions).astype(float)
        except:
            error('Basis positions must be an N element list of 3D coordinates')

        self.num_basis = self.basis_positions.shape[0]
        if self.basis_positions.shape[1] != 3:
            error('Must provide 3D coordinates for all basis atoms')
        self.reduced_coords = bool(reduced_coords)

        self.basis_types = np.array(basis_types).astype(str)
        if self.basis_types.shape[0] != self.basis_positions.shape[0]:
            error('Number of atom types must match the number of positions given')
        self.unique_types = list(np.unique(self.basis_types))
        self.num_types = len(self.unique_types)

        if masses != False:
            try:
                self.masses = np.array(masses).astype(float)
            except:
                error('Masses must be floats')
            if len(np.argwhere(self.masses < 0.0001)) != 0:
                error('Masses must be positive')
            if len(self.masses) != self.num_types:
                error('Must provide 1 mass per atom type')
        else:
            self.masses = np.ones(self.num_types)

        self.basis_index = np.arange(1,self.num_basis+1).reshape(self.num_basis,1)
        self.unitcell_index = np.ones((self.num_basis,1))
        self.mass_list = np.zeros((self.num_basis,1))
        for i in range(self.num_types):
            self.mass_list[np.argwhere(self.basis_types == 
                self.unique_types[i])] = self.masses[i]

    def lattice_vectors(self,lattice_vectors,lattice_constants=[1,1,1]):
        try:
            self.lattice_vectors = np.array(lattice_vectors).astype(float)
        except:
            error('Lattice vectors must be a 3 element list of 3D vectors')
        if self.lattice_vectors.shape[0] != 3:
            error('Expected 3 lattice vectors')
        if self.lattice_vectors.shape[1] != 3:
            error('Lattice vectors must be 3D')
        try:
            self.lattice_constants = np.array(lattice_constants).astype(float)
        except:
            error('Lattice constants must be a 3 element list of floats')

        self.lattice_vectors[0,:] = self.lattice_vectors[0,:]*self.lattice_constants[0]
        self.lattice_vectors[1,:] = self.lattice_vectors[1,:]*self.lattice_constants[1]
        self.lattice_vectors[2,:] = self.lattice_vectors[2,:]*self.lattice_constants[2]
        if self.reduced_coords == True:
            for i in range(self.basis_positions.shape[0]):
                self.basis_positions[i,:] = (
                        self.basis_positions[i,0]*self.lattice_vectors[0,:]+
                        self.basis_positions[i,1]*self.lattice_vectors[1,:]+
                        self.basis_positions[i,2]*self.lattice_vectors[2,:])

    def replicate(self,num_reps=[1,1,1]):
        self.num_reps = num_reps
        atom_types = np.copy(self.basis_types)
        pos = np.copy(self.basis_positions[:,:])
        basis = np.copy(self.basis_index)
        mass = np.copy(self.mass_list)
        unitcells = np.copy(self.unitcell_index)

        for i in range(num_reps[0]):
            a1 = np.copy(self.lattice_vectors[0,:])*i
            for j in range(num_reps[1]):
                a2 = np.copy(self.lattice_vectors[1,:])*j
                for k in range(num_reps[2]):
                    a3 = np.copy(self.lattice_vectors[2,:])*k
                    if i == 0 and j == 0 and k == 0: 
                        continue
                    new_cell = np.copy(self.basis_positions)
                    new_cell[:,0] = new_cell[:,0]+a1[0]+a2[0]+a3[0]
                    new_cell[:,1] = new_cell[:,1]+a1[1]+a2[1]+a3[1]
                    new_cell[:,2] = new_cell[:,2]+a1[2]+a2[2]+a3[2]

                    unitcells = unitcells+1
                    self.unitcell_index = np.append(self.unitcell_index,unitcells,axis=0)
                    self.basis_index = np.append(self.basis_index,basis,axis=0)
                    self.mass_list = np.append(self.mass_list,mass,axis=0)

                    pos = np.append(pos,new_cell,axis=0)
                    atom_types = np.append(atom_types,self.basis_types,axis=0)

        pos[:,0] = pos[:,0]-pos[:,0].min()
        pos[:,1] = pos[:,1]-pos[:,1].min()
        pos[:,2] = pos[:,2]-pos[:,2].min()
        self.atom_positions = pos
        self.atom_types = atom_types
        self.num_atoms = pos.shape[0]

    def write_abipos(self,filename='structure.abi'):
        self.unique_types = list(np.unique(self.atom_types))
        self.num_types = len(self.unique_types)
        with open(filename,'w') as fid:
            fid.write(' # ')
            for i in range(self.num_types):
                fid.write('{:d} = {}, '.format(i+1,self.unique_types[i]))
            fid.write('\n\n typat   ')
            for i in range(self.num_atoms):
                fid.write('{:d} '.format(self.unique_types.index(self.atom_types[i])+1))
            fid.write('\n\n xangst   {:10f} {:10f} {:10f}\n'
                    .format(self.atom_positions[0,0],self.atom_positions[0,1],
                        self.atom_positions[0,2]))
            for i in range(1,self.num_atoms):        
                fid.write('          {:10f} {:10f} {:10f}\n'
                        .format(self.atom_positions[i,0],self.atom_positions[i,1],
                        self.atom_positions[i,2]))

    def write_xyz(self,filename='structure.xyz'):
        with open(filename,'w') as fid:
            fid.write('{:d}\n{}\n'.format(self.num_atoms,self.description))
            for i in range(self.num_atoms-1):
                fid.write('{} {:.10f} {:.10f} {:.10f}\n'
                        .format(self.atom_types[i],
                                self.atom_positions[i,0],
                                self.atom_positions[i,1],
                                self.atom_positions[i,2]))
            fid.write('{} {:.10f} {:.10f} {:.10f}\n'
                        .format(self.atom_types[-1],
                                self.atom_positions[-1,0],
                                self.atom_positions[-1,1],
                                self.atom_positions[-1,2]))

    def write_lammps(self,filename='structure.lammps'):
        if (self.lattice_vectors[0,1] != 0 or self.lattice_vectors[0,2] != 0 or
                self.lattice_vectors[1,0] != 0 or self.lattice_vectors[1,2] != 0 or
                self.lattice_vectors[2,0] != 0 or self.lattice_vectors[2,1] != 0):
            error('Can only create LAMMPS file for orthogonal lattice vectors')
            exit()

        with open(filename,'w') as fid:
            fid.write(self.description+'\n\n')
            fid.write('{} atoms\n\n'.format(self.num_atoms))
            fid.write('{} atom types\n\n'.format(self.num_types))
            fid.write('{:12.8f} {:12.8f} xlo xhi\n'.format(0,
                self.lattice_vectors[0,0]*self.num_reps[0]))
            fid.write('{:12.8f} {:12.8f} ylo yhi\n'.format(0,
                self.lattice_vectors[1,1]*self.num_reps[1]))
            fid.write('{:12.8f} {:12.8f} zlo zhi\n\nMasses\n\n'
                    .format(0,self.lattice_vectors[2,2]*self.num_reps[2]))
            for i in range(self.num_types):
                fid.write('{} {:.4f}\n'.format(i+1,self.masses[i]))
            fid.write('\nAtoms\n\n')
            for i in range(self.num_atoms-1):
                fid.write('{} {} {:12.8f} {:12.8f} {:12.8f}\n'
                        .format(i+1,self.unique_types.index(self.atom_types[i])+1,
                        self.atom_positions[i,0],self.atom_positions[i,1],
                        self.atom_positions[i,2]))
            fid.write('{} {} {:12.8f} {:12.8f} {:12.8f}'
                        .format(i+2,self.unique_types.index(self.atom_types[-1])+1,
                        self.atom_positions[-1,0],self.atom_positions[-1,1],
                        self.atom_positions[-1,2]))

    def write_lattice_file(self,filename='lattice.dat'):
        data = np.append(np.arange(1,self.num_atoms+1).reshape(self.num_atoms,1),
                self.unitcell_index,axis=1)
        data = np.append(data,self.basis_index,axis=1)
        data = np.append(data,self.mass_list,axis=1)
        np.savetxt(filename,data,fmt='%d %d %d %.4f')






