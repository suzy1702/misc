import numpy as np

def write_output(params,data):
    np.savetxt(params.out_prefix+'.EIGS',data.eigenvals) #,fmt='%4f')
def read_eigs(params):
    eigenvals = np.loadtxt(params.out_prefix+'.EIGS')
    return eigenvals
