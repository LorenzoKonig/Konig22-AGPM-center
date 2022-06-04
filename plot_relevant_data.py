# Plot the results of the MEEP simulations. By Lorenzo König 2022.
# =========================================================================

# This script plots the three relevant plots, reading the relevant data
# from 'relevant_data.npy.npz', previously saved in 'vortex_mask.py'. The 
# data was saved in meep without pml layers. 

# -------------------------------------------------------------------------
# general parameters

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
from matplotlib.patches import Rectangle

design = 'AGPM'      #  specify design name
paths  = ['.']       #  paths where to look for simulation results
wvls   = ['3.5000']  #  plot results for these wavelengths
nulls  = []          #  initialize null dpeth list

res = 20    #  specify relevant simulation parameters - adjust plot limits
gp  = 1.21  #  accordingly

rel     = np.load('relevant_data.npy.npz')
ee      = rel['arr_0']
eer     = rel['arr_1']
eernorm = rel['arr_2']  #  transmitted intensity (LHC) - use this to normalize leakage
eel     = rel['arr_3']

eer = eer/(np.mean(eernorm)+np.mean(eer))  #  normalize leakage

for path in paths:
  for wvl in wvls:
    vminrl=-4  #  colorscale limits
    vmaxrl=0   #
    
    fig=plt.figure(figsize=[14.,6.])
    
    plt.subplot(131)  #  structure
    cmm=plt.get_cmap('binary',2)
    im_eps=plt.imshow(ee.transpose(),cmap=cmm,vmin=1,vmax=2.38**2)
    plt.gca().invert_yaxis()
    plt.ylabel("y [\u03BCm]")  #  add y-labels
    plt.xlabel("x [\u03BCm]")  #  add x-labels
    plt.xticks([0,10*res,20*res,30*res],[0,10,20,30])
    plt.yticks([0,10*res,20*res,30*res],[0,10,20,30])

    plt.subplot(132)  #  leakage
    im_int=plt.imshow(np.log10(eer).transpose(),cmap='hot',vmin=vminrl,vmax=vmaxrl)
    plt.gca().invert_yaxis()
    plt.xlabel("x [\u03BCm]")  #  add x-labels
    plt.xticks([0,10*res,20*res,30*res],[0,10,20,30])
    plt.yticks([0,10*res,20*res,30*res],[0,10,20,30])
    
    plt.subplot(133)  #  phase ramp
    im_phi=plt.imshow(eel.transpose(),cmap='hsv',vmin=-np.pi,vmax=np.pi)
    plt.gca().invert_yaxis()
    plt.xlabel("x [\u03BCm]")  #  add x-labels
    plt.xticks([0,10*res,20*res,30*res],[0,10,20,30])
    plt.yticks([0,10*res,20*res,30*res],[0,10,20,30])
    
    cbar_ax=fig.add_axes([.125,.8,.228,.020])     # colorbar for design plots
    cl_eps=fig.colorbar(im_eps,cax=cbar_ax,ticks=[2.1661,4.4983],orientation='horizontal')
    cl_eps.set_label('Refractive index n',labelpad=10)
    cl_eps.ax.set_xticklabels(['1.00','2.38'])
    cl_eps.ax.xaxis.set_ticks_position('top')
    cl_eps.ax.xaxis.set_tick_params(length=0,width=0,which='major')
    cl_eps.ax.xaxis.set_label_position('top')

    cbar_ax=fig.add_axes([.399,.8,.228*4.706/4.,.020]) # colorbar for intensity plots # max is at 10^-.600, resulting in 4.706/4 factor
    cl_int=fig.colorbar(im_int,cax=cbar_ax,ticks=[-4,-3,-2,-1],orientation='horizontal')  #
    cl_int.set_label('Leakage          ',labelpad=10)
    cl_int.ax.set_xticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$'])
    cl_int.ax.xaxis.set_ticks_position('top')
    cl_int.ax.xaxis.set_label_position('top')
    
    cbar_ax=fig.add_axes([.6725,.8,.228,.020])     # colorbar for phase plots
    cl_phi=fig.colorbar(im_phi,cax=cbar_ax,ticks=[-np.pi,0,np.pi],orientation='horizontal')
    cl_phi.set_label('Phase',labelpad=10)
    cl_phi.ax.set_xticklabels(['-$\pi$','0',"$\pi$"])
    cl_phi.ax.xaxis.set_ticks_position('top')
    cl_phi.ax.xaxis.set_label_position('top')

    # to make the colorbar end at the max value uncoment the following lines and change the relevant colorbar properties (ticks etc.)
    axxx=fig.add_axes([.399+.228,.8-.020*.155,.228*.706/4.*1.05,.020*1.31])
    axxx.add_patch(Rectangle((0,0),1,1,facecolor='white',edgecolor='white'))
    axxx.axis('off')
    axxx.axvline(x=0.,ymin=0.15,ymax=.85,color='k',lw=1.5)
    
    # to make design colorbar as two small rectangles uncomment the following
    axxx1=fig.add_axes([.125-.02,.8-.020*.155,.228/4-.004+.02,.020*1.31])
    axxx1.add_patch(Rectangle((0,0),1,1,facecolor='white',edgecolor='white'))
    axxx1.axis('off')
    axxx1.axvline(x=1.,ymin=0.15,ymax=.85,color='k',lw=1.5)
    axxx2=fig.add_axes([.125+.228*1/4+.004,.8-.020*.155,.228*1/2-2*.004,.020*1.31])
    axxx2.add_patch(Rectangle((0,0),1,1,facecolor='white',edgecolor='white'))
    axxx2.axis('off')
    axxx2.axvline(x=0.,ymin=0.15,ymax=.85,color='k',lw=1.5)
    axxx2.axvline(x=1.,ymin=0.15,ymax=.85,color='k',lw=1.5)
    axxx3=fig.add_axes([.125+.228*3/4+.004,.8-.020*.155,.228*1/4,.020*1.31])
    axxx3.add_patch(Rectangle((0,0),1,1,facecolor='white',edgecolor='white'))
    axxx3.axis('off')
    axxx3.axvline(x=0.,ymin=0.15,ymax=.85,color='k',lw=1.5)

    # to save as png uncomment the following
    #plt.savefig('results_vortex_mask_'+design+'_'+wvl+'.png',bbox_inches='tight',pad_inches=.2)

    # to save as pdf uncomment the following
    plt.savefig('results_vortex_mask_'+design+'_'+wvl+'.pdf',bbox_inches='tight',pad_inches=.2,dpi=190) # set the dpi to >~120 to solve a known bug in plt.savefig when saving to pdf (this made the colorbar not to be centred in its axes frame); https://github.com/matplotlib/matplotlib/issues/6827

    plt.close('all')
