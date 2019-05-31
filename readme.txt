Matlab files to generate figure 6 from O'Donnell et al, J Neurosci, 2011.

Code simulates 3-compartment model: dendritic spine head, dendritic
spine neck, and adjoining dendritic branch. Simulation includes both
electrical and Ca2+ dynamics. See O'Donnell et al paper Methods
section and Figure 6 for further details.

Included are two files:
- spinemodel_fig6.m specifies all model and simulation parameters, calls
  spineodes.m, solves differential equations using a Matlab ODE solver,
  and plots figures.
- spineodes.m specifies system ODEs for use by spinemodel_fig6.m.

To try a basic simulation, simply run spinemodel_fig6.m as is from Matlab.

To vary model parameters, simulation parameters or stimulation
protocol, open spinemodel_fig6.m and investigate. Parameters are
clearly notated.

Any queries or comments can be sent to cian@salk.edu.
