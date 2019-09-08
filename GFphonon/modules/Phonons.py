"""
TO DO:

2) include therodynamic beta in the Greens functions

3) add acoustic sum rule to FC matrix

4) currently just simple time averaging the displacements. should consider how to
    calculate the ensemble averaging correctly; if its corr. fx, best way might
    be IFFT(FFT(a)*FFT(b).conj())

5) average across 'ensembles'

"""

import numpy as np
import scipy.linalg as linalg

class phonons:
    def compute_greens_functions(self,params,data):
        if params.debug:
            num_loops = 1
        else:
            num_loops = params.num_ensembles
        for i in range(num_loops):
            self.loop_index = i
            print('\nNow on ensemble {} out of {}\n'.format(i+1,num_loops))
            self.get_simulation_data(params,data)
            self.compute_normal_mode(params,data)
            data.gf_total = data.gf_total+data.gf
        data.gf_total = data.gf_total/num_loops 
        self.diagonalize(params,data)

    def compute_normal_mode(self,params,data):
        for q in range(params.num_qpoints):
            print('\tq = {} out of {}'.format(q+1,params.num_qpoints))
            self.q_index = q
            for b in range(params.num_basis):
                self.basis_index = b
                self.qpoint = data.qpoints[q,:]
                self.compute_qfac(params,data)
                self.pos = data.pos[:,np.argwhere(params.basis_pos == (b+1)).flatten(),:]
                self.pos = (self.pos*self.qfac).sum(axis=1)/np.sqrt(params.num_unitcells)
                data.qpos[:,b,:] = self.pos
            self.compute_matrix_elements(params,data)

    def compute_qfac(self,params,data):
        self.qfac = np.zeros((params.num_unitcells,1)).astype(complex)
        self.qfac = np.exp(-1j*(data.cell_vecs*
            np.tile(self.qpoint,(params.num_unitcells,1)))).sum(axis=1)
        self.qfac = np.tile(self.qfac.reshape(params.num_unitcells,1),
                (params.steps_per_ensemble,1,3))     

    def compute_matrix_elements(self,params,data):
        for k in range(params.num_basis):
            xmass = params.masses[k]
            for a in range(3):
                for l in range(params.num_basis):
                    ymass = params.masses[l]
                    for b in range(3):
                        data.gf[self.q_index,k*3+a,l*3+b] = (
                            (data.qpos[:,k,a]*data.qpos[:,l,b].conj()).mean()-
                            (data.qpos[:,k,a].mean()*data.qpos[:,l,b].conj().mean())
                            )/np.sqrt(xmass*ymass)
    
    def diagonalize(self,params,data):
        for q in range(params.num_qpoints):
            data.eigenvals[q,:], vecs = linalg.eigh(linalg.inv(data.gf_total[q,:,:]))

    def get_simulation_data(self,params,data):
        data.pos = data.database['pos'][self.loop_index*params.steps_per_ensemble:
                (self.loop_index+1)*params.steps_per_ensemble,:,:]
        data.cell_vecs = data.pos[:,np.argwhere(params.basis_pos == 1)
                .flatten(),:].mean(axis=0)
