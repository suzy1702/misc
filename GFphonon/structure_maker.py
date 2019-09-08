import sys
import numpy as np
modulepath = './modules'
sys.path.append(modulepath)
import Lattice


crystal = Lattice.structure_maker('Si crystal')

crystal.basis(
positions = [[0.0, 0.0, 0.0],
             [0.0, 1/2, 1/2],
             [1/2, 0.0, 1/2],
             [1/2, 1/2, 0.0],
             [1/4, 1/4, 1/4],
             [1/4, 3/4, 3/4],
             [3/4, 1/4, 3/4],
             [3/4, 3/4, 1/4]],
basis_types = ['1',
               '1',
               '1',
               '1',
               '1',
               '1',
               '1',
               '1'],
masses =   [28.0855],
reduced_coords=True)

crystal.lattice_vectors(
lattice_vectors =   [[ 1, 0, 0],
                     [ 0, 1, 0],
                     [ 0, 0, 1]],
lattice_constants = [5.431, 5.431, 5.431])

crystal.replicate([16,2,2])

crystal.write_xyz('structure.xyz')
crystal.write_lammps('data.lammps')
crystal.write_lattice_file('lattice.dat')

