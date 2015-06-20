# GENSO
A tool for generating or refining kinetic descriptions of seismic rupture plane

This software take an existing slip distribution and refines the descritisation. Those kinematic descriptions are used for simulating ground motions from earthquakes. The earthquake source is descritised into a number of patches, for each of which the position, the slip, the fault mechanism, the rupture time and the rise time are specified.

Those kinematic source descriptions are usually found from inversion of seismic are geodetic data. Unfortunately the resolution of the inversions is usually poor. For the forward simulation of seismograms, much finer discretisations are needed to avoid spatial aliasing.

GENSO reads an exisiting source description, refines the discretisation to the specified patch-size and adds heterogeneities in slip and fault mechanism such that the final distributions are self-similar. The characteristics of the distributions at low wavenumbers are kept.
 
- generates complete source description for a kinematic source model
- reads input from file "infile" (user defined) 
- outputs one file XXX_genso for each source segment input which includes refined slip distribution
- outputs the final source description in the file source.out including the refines slip, 
       patch rupture times, patch rise times

- comment: rise time is defined as 1/(2 pi fc) , where fc is the corner frequency of the source time function.
           The rise time is the time from rupture initiation to the peak moment rate.
           this approach makes use of Brune's source time function
  
- comment2: The rise time is scaled such that the radiated seismic energy of the simulated rupture matches the energy magnitude 
                  specified in the input file
