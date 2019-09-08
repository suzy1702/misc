import numpy as np
import matplotlib.pyplot as plt
import FileIO

def plot_bands(params):
    eigenvals = FileIO.read_eigs(params)

    fig, ax = plt.subplots()
    fig.set_size_inches(5,6,forward=True)
    fig.tight_layout(pad=5)
    
    for v in range(params.num_basis*3):
        ax.plot(eigenvals[:,v],ls='',lw=1,ms=5,mfc='b',mec='k',mew=1,marker='o')

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)

    ax.minorticks_on()
    ax.tick_params(which='both', axis='y', width=1, labelsize='large')
    ax.tick_params(which='both', axis='x', width=1, labelsize='small',
            labelrotation=0.0,pad=5.0)
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=3, color='k')
#    ax.tick_params(axis='x',which='both',labelbottom=False)
#    ax.set_xlabel(r'$\bfq$',labelpad=5.0,fontweight='normal',fontsize='x-large')
#    ax.set_ylabel(r'$\omega$ (THz)',labelpad=3.0,fontweight='normal',fontsize='x-large')
#    fig.suptitle(r'$\Phi$($\bfq$,$\omega)$',y=0.95,fontsize='x-large')
#    plt.savefig('example.png',format='png',dpi=300,bbox_inches='tight')
    plt.show()

