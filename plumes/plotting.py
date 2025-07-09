import numpy as np
import sys,pdb
from scipy.special import erfcinv as erfcinv
from calc_sigmas import calc_sigmas 
import matplotlib.pyplot as plt

def mole_fraction_plot(mole_fraction=None,x=None,y=None,x_slices=[0],title='CH4',plot_lims=[2,3],cmap=plt.cm.hot_r,xlabel='Downwind Direction (km)',ylabel='Crosswind Direction (km)'):
    fig,axs = plt.subplots(1,2,figsize=(12,4))
    ax = axs[0]
    g=ax.pcolormesh(x/1000,y/1000,mole_fraction[:,:],cmap=cmap,vmin=plot_lims[0],vmax=plot_lims[1]); plt.colorbar(g,ax=ax)
    for xs in x_slices:
        ind = np.argmin((x-xs)**2)
        ax.plot(np.array([xs/1000,xs/1000]),np.array([y.min()/1000,y.max()/1000]),'--');
    ax.set_ylabel(ylabel,fontsize=16)
    ax.set_xlabel(xlabel,fontsize=16)
    ax.set_title(f'{title}',fontsize=18)

    ax = axs[1]
    g = []
    l = []
    for xs in x_slices:
        ind = np.argmin((x-xs)**2)
        g0=ax.plot(mole_fraction[:,ind],y/1000);
        g.append(g0[0])
        l.append(f'{xs/1000}km')

    ax.legend(g,l)
    ax.set_ylabel(ylabel,fontsize=16)
    ax.set_xlabel('Mole Fraction (ppm)',fontsize=16)
    ax.set_xlim(plot_lims)
    ax.set_title(f'Slice Profile of \n {title}',fontsize=18)
    return
