# GENSO
A tool for generating or refining kinetic descriptions of seismic rupture plane
 
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
