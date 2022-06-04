# Konig22-AGPM-center
Supporting material for König et al., "Optimal Design of the Annular Groove Phase Mask Central Region" (2022). 

The code presented in this repository produces the simulation results for the AGPM center (Figure 5 in the paper). Other grating geometries, such as straight grating lines (Figure 3 in the paper), are straightforward to implement. 

To execute the code you will have to install MEEP (https://meep.readthedocs.io/en/latest/) using the python interface. 

Use 'vortex_mask.py' to run the MEEP simulation passing the simulation parameters via command line arguments (python vortex_mask.py -resolution=20 -dpml=12 etc.). 
Use 'plot_relevant_results.py' to plot the design, the leakage, and the phase ramp. 
