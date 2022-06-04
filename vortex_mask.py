# MEEP/FDTD simulation of Annular Groove Phase Mask. By Lorenzo König 2022.
# =========================================================================

# Meep is a free and open-source software package for electromagnetics 
# simulation via the finite-difference time-domain (FDTD) method. More
# information about MEEP can be found in [1] and on 
# https://meep.readthedocs.io/en/latest/ . 
# [1] A. Oskooi, D. Roundy, M. Ibanescu, P. Bermel, J.D. Joannopoulos, and
#     S.G. Johnson, “MEEP: A flexible free-software package for 
#     electromagnetic simulations by the FDTD method,” Computer Physics 
#     Communications, Vol. 181, pp. 687-702 (2010).

# MEEP simulation of AGPM
# =======================

# This script propagates a circularly polarized plane wave through a Vortex
# Phase Mask consiting of full-dielectric (diamond) subwavelength 
# gratings. The simulation cell is composed of (top-bottom): air padding - 
# source - air padding - diamond MS structure - diamond substrate - 
# air padding . The electric fields are monitored below the substrate (and 
# possibly before and inside the MS structure, inside the substrate, and at 
# different distances below the substrate) to retrieve the phaseramp and 
# leakage. UPDATE: the air padding on the substrate side has been filled 
# with diamond, because only the grating behaviour is relevant. 

import meep as mp
import numpy as np
import argparse
import matplotlib.pyplot as plt

print("MEEP simulation of AGPM\n" \
      "---\n" \
      "parameters to be set:\n" \
      " -resolution=..\n" \
      " -wvl=..        (wavelength of monochromatic source)\n" \
      " -sourc=..      (source; 4 is default RHC input)\n" \
      " -turnontime=.. (time to turn the source on)\n" \
      " -dpml=..       (pml thickness)\n" \
      " -dpad=..       (air padding above structure)\n" \
      " -dsub=..       (substrate thickness)\n" \
      " -dback=..      (air padding below substrate)\n" \
      " -gp=..         (grating period of AGPM)\n" \
      " -gh=..         (height of grating/MS structure)\n" \
      " -gdc=..        (grating duty cycle of AGPM)\n" \
      " -dpillar=..    (diameter of central pillar of AGPM; in lw)\n" \
      " -num_cells=..  (# of grating lines of AGPM; x2 for diameter)\n" \
      " -rununtil=..   (runtime of simulation)\n" \
      " -design=..     (choose MS structure design name)\n" \
      " -twos=..     (output only two slices? [True/False])\n")

# -------------------------------------------------------------------------
# parse arguments from command line

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-resolution', type=float, default=10, 
     help='resolution (default: 10 pxl/μm')
  parser.add_argument('-wvl', type=float, default=3.5, 
     help='source wavelength (default: 3.5 = center of METIS L band)')
  parser.add_argument('-sourc', type=int, default=4, 
     help='sourc (default: 4 = RHC)')
  parser.add_argument('-turnontime', type=float, default=0.0, 
     help='turnon time for the source to reduce oscillations (default: 0)')
  parser.add_argument('-dpml', type=float, default=4.0, 
     help='pml thickness (default: 4.0 μm)')
  parser.add_argument('-dpad', type=float, default=3.0, 
     help='padding (air thickness) (default: 3.0 μm)')
  parser.add_argument('-dsub', type=float, default=3.0, 
     help='substrate thickness (default: 3.0 μm)')
  parser.add_argument('-dback', type=float, default=2.0, 
     help='air thickness below AGPM (default: 2.0 μm)')
  parser.add_argument('-gp', type=float, default=1.21, 
     help='grating periodicity (default: 1.21 μm)')
  parser.add_argument('-gh', type=float, default=4.5753, 
     help='grating height (default: 4.5753 μm)')
  parser.add_argument('-swa', type=float, default=0., 
     help='sidewall angle in deg (default: 0 deg)')
  parser.add_argument('-gdc', type=float, default=0.6509, 
     help='grating duty cycle (default: 0.6509)')
  parser.add_argument('-dpillar', type=float, default=2.0, 
     help='diameter of central pillar in units of linewidth (default: 2)')
  parser.add_argument('-num_cells', type=int, default=10, 
     help='number of grating lines (default: 10)')
  parser.add_argument('-rununtil', type=int, default=100, 
     help='run until some time (default: 100)')
  parser.add_argument('-design', type=str, default='AGPM', 
     help='name of design (default: AGPM)')
  parser.add_argument('-twos', type=bool, default=False, 
     help='output only two zslices? (Default: False)')
  args = parser.parse_args()
  print(args)

# -------------------------------------------------------------------------
# define the simulation cell

#  /_________________ /    #  UPDATE: don't use air layer after the AGPM, 
# |  ______pml_____  |     #  rather fill the 'dback' layer and the pml on
# |.|.....s.r.c....|.|     #  that side with diamond (what matters is the
# |_|_____A_I_R____|_|     #  effect of the grating, not the back side). 
# |_|XXX_VORTEX_XXX|_|     #
# | |     S U B    | |     #    z  y
# |.|.....m.o.n....|.|     #         
# | |_____S_U_B____| |     #    | /  
# |________pml_______|/    #    |/__ x

design = args.design  #  vortex mask design

resolution = args.resolution  #  pxl/um - 10 should be fine
dpml       = args.dpml        #  PML thickness - 8 is good (4 ok for tests)
dpad       = args.dpad        #  air thickness above vortex mask
gh         = args.gh          #  height of MS strucutre
dsub       = args.dsub        #  substrate thickness
dback      = args.dback       #  air thickness below substrate

gp = args.gp                      # grating period
gdc = args.gdc                    # filling factor (grating duty cycle)
swa = args.swa*np.pi/180          # side wall angle
num_cells = args.num_cells        # number of grating lines to simulate
sx = dpml+(num_cells*gp*2)+dpml   # re-define cell size in case of AGPM
sy = sx                           #

sz = dpml+dpad+gh+dsub+dback+dpml  #  cell size

cell_size = mp.Vector3(sx,sy,sz)  #  define cell size

pml_layers = [mp.PML(thickness=dpml)]  # make pml's on all sides

# -------------------------------------------------------------------------
# define the source

wvl = args.wvl  #  monochromatic source wavelength

sourc = args.sourc  # 'source' is builtin in python, so use 'sourc'
if sourc==3:
    a_x = 1
    a_y = 0-1j
elif sourc==4:
    a_x = 1
    a_y = 0+1j
else:
    print('sourc must be 3 (LHC pw) or 4 (RHC pw)! Stopping script.')
    exit() # stop if sourc is not 3 or 4 (circularly polarized plane wave)

turnontime = args.turnontime

sources = [mp.Source(mp.ContinuousSource(frequency=1/wvl,width=turnontime),
                component=mp.Ex,
                center=mp.Vector3(0,0,0.5*sz-dpml-0.5*dpad),
                size=mp.Vector3(sx,sy,0),
                amplitude=a_x),
           mp.Source(mp.ContinuousSource(frequency=1/wvl,width=turnontime),
                component=mp.Ey,
                center=mp.Vector3(0,0,0.5*sz-dpml-0.5*dpad),
                size=mp.Vector3(sx,sy,0),
                amplitude=a_y)]

# -------------------------------------------------------------------------
# define the geometry

air = mp.Medium(index=1.00)      #  index of refraction
diamond = mp.Medium(index=2.38)  #

geometry = [mp.Block(material=diamond,  #  substrate block
                     size=mp.Vector3(mp.inf,mp.inf,dsub+dback+dpml),
                     center=mp.Vector3(0,0,-0.5*sz+0.5*(dpml+dback+dsub)))]

if design=='AGPM':
    # AGPM design consists of concentric grating lines
    dpillar = args.dpillar  #  diameter of central pillar
    jmax = np.int(np.floor(num_cells*np.sqrt(2)+2)+np.floor(np.sqrt(2)*dpml/gp))
    if dpillar<=0:
        jmax = jmax-1 # use -1 in the following for loop to avoid ..
        # .. drawing the last diamond cylinder which may have negative ..
        # .. radius; draw the last air cylinder separately outside for loop
    rmax = np.floor((sx/2-dpml)/gp*np.sqrt(2)+2)+np.floor(np.sqrt(2)*dpml/gp) # max radius of grating ..
    # .. lines to be used in geometry (must start from outermost lines ..
    # .. since solid cylinder objects alternating air and diamond are used)
    for j in range(jmax):
        geometry.append(mp.Cone(material=air,
                     center=mp.Vector3(0,0,sz/2-dpml-dpad-gh/2),
                     radius=rmax*gp-j*gp-(2-dpillar)/2*gdc*gp-gh*np.tan(swa),
                     radius2=rmax*gp-j*gp-(2-dpillar)/2*gdc*gp,
                     height=gh))
        geometry.append(mp.Cone(material=diamond,
                     center=mp.Vector3(0,0,sz/2-dpml-dpad-gh/2),
                     radius=rmax*gp-j*gp-(1-gdc)*gp-(2-dpillar)/2*gdc*gp+gh*np.tan(swa),
                     radius2=rmax*gp-j*gp-(1-gdc)*gp-(2-dpillar)/2*gdc*gp,
                     height=gh))
    if dpillar<=0:
        geometry.append(mp.Cone(material=air, # jmax was already decreased by 1, so no -1 here
                     radius=rmax*gp-jmax*gp-(2-dpillar)/2*gdc*gp-gh*np.tan(swa),
                     radius2=rmax*gp-jmax*gp-(2-dpillar)/2*gdc*gp,
                     height=gh))

else:
    print('Design not supported. Stopping script.')
    exit() # stop if design is not AGPM

# other deisngs can be implemented by defining the appropriate geometry

# -------------------------------------------------------------------------
# define simulation

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    geometry=geometry,
                    sources=sources,
                    force_complex_fields=True)

# -------------------------------------------------------------------------
# run simulation and get results

rununtil = args.rununtil  #  simulation runtime, in meep units

twos = args.twos  #  output only two slices?
if twos==True:  # output only one slice inside substrate and one below
   zslices = [-sz/2+dpml+dback+dsub*1/2,-sz/2+dpml+dback/2]
else:
   zslices = [sz/2-dpml-dpad*2/3,      #  output E-fields in these z-planes
           sz/2-dpml-dpad-gh/2,        #
           -sz/2+dpml+dback+dsub*1/2,  #
           -sz/2+dpml+dback*3/4,       #
           -sz/2+dpml+dback*2/4,       #
           -sz/2+dpml+dback*1/4]       #

eps_data = len(zslices)*[0]
ex_data = len(zslices)*[0]
ey_data = len(zslices)*[0]
def get_zslices(sim):
    for i in range(len(zslices)):
        eps_data[i] = sim.get_array(center=mp.Vector3(0,0,zslices[i]),
                                    size=mp.Vector3(sx-2*dpml,sy-2*dpml,0),
                                    component=mp.Dielectric)
        ex_data[i] = sim.get_array(center=mp.Vector3(0,0,zslices[i]),
                                    size=mp.Vector3(sx-2*dpml,sy-2*dpml,0),
                                    component=mp.Ex)
        ey_data[i] = sim.get_array(center=mp.Vector3(0,0,zslices[i]),
                                    size=mp.Vector3(sx-2*dpml,sy-2*dpml,0),
                                    component=mp.Ey)

sim.run(mp.at_end(get_zslices),until=rununtil)

# -------------------------------------------------------------------------
# save all slices

np.save('eps_'+design+'_{:.4f}'.format(wvl)+'.npy',eps_data)  #  1st dimension is z, then x and y
np.save('ex_'+design+'_{:.4f}'.format(wvl)+'.npy',ex_data)    #
np.save('ey_'+design+'_{:.4f}'.format(wvl)+'.npy',ey_data)    #

np.save('params_'+design+'_{:.4f}'.format(wvl)+'.npy',vars(args))  #  saves a dictionary to a np.array
#  (  read from parameter file like this: 
#                   data_array = np.load('parameters.npy',allow_pickle=True)
#                   dictionary = data_array.flatten()[0]
#                   resolution = dictionary['resolution']                  )

# -------------------------------------------------------------------------
# save only relevant results

eps_data = np.moveaxis(eps_data,0,-1) #  zxy
ex_data = np.moveaxis(ex_data,0,-1)   #  -->
ey_data = np.moveaxis(ey_data,0,-1)   #  xyz

er_data=1./np.sqrt(2)*(ex_data-1j*ey_data)
el_data=1./np.sqrt(2)*(ex_data+1j*ey_data)

er_int=np.abs(er_data)**2  #  calculate phase and intensity of R/L
er_phi=np.angle(er_data)   #  polarizations
el_int=np.abs(el_data)**2  #
el_phi=np.angle(el_data)   #

if twos==True:
  ee = eps_data[:,:,0]
else:
  ee = eps_data[:,:,1]
eer = er_int[:,:,-1]
eernorm = el_int[:,:,-1]
eel = el_phi[:,:,-1]

np.savez('relevant_data.npy',ee,eer,eernorm,eel)

