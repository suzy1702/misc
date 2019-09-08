import Lattice
import h5py
import numpy as np

class data(Lattice.lattice):
    def structures(self,params):
        data.eigenvals = np.zeros((params.num_qpoints,params.num_basis*3))
        self.database = h5py.File(params.database_file,'r')
        self.qpos = np.zeros((params.steps_per_ensemble,
            params.num_basis,3)).astype(complex)
        self.gf = np.zeros((params.num_qpoints,params.num_basis*3,
            params.num_basis*3)).astype(complex)
        self.gf_total = np.zeros((params.num_qpoints,params.num_basis*3,
            params.num_basis*3)).astype(complex)


