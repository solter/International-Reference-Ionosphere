# SOURCE F90 FILES
This directory contains f90 source files.

They begin by running the original f77 files
through the refactoring scripts (fixLineCont.awk
then convert2indent.awk). 
At this point they were renamed to name.0.f90

## Refactoring state
The number between the root name and the f90
extension represent the state of the refactor.
Zero indicates that no manual steps at refactoring have
been performed, only the awk scripts.

Manual refactoring steps indicate the following:

1. Manual updates to remove all gotos

## Files
The are as follows:
* *cira* - Neutral temperature and density model
* *igrf* - calculates international geomagnetic reference field
* *iridreg* - d region model by Friedrich and Torkar
* *iriflip* - calculates ion densities via the FLIP (IDC?) model
* *irifun* - a hodgepodge of functions iri uses. 
  This has been broken into the following:
  * *elec_dens* - electron density functions
  * *epstein_lay* - Epstein and Lay functions
  * *ion_dens* - ion density functions
  * *peaks* - peak parameters (HmF2, FoF2, NmE, etc) functions
  * *profile* - ionospheric profile functions
  * *temperature* - temperature calculations
  * *time* - time calculations
  * *vdrift* - vertical drift calculations
  * *irifun* - functions that did not fit in the above categories
* *irisub* - the main routines to interface with the IRI model.
             essentially a library containing 2 subroutins
  * irisub - calculates IRI profile at 1 lat/lon
  * iriweb - cycles through a parameter set calculating the profiles
* *iritec* - routines to calculate the total electron content
* *iritest* - a sample test routine. Inputs and outputs are via interactive user inputs.
