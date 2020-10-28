***************************
Try to make DNA homogeneous
***************************

The data from :expt:`46` and :expt:`68` is make harder to interpret because the 
DNA is not particularly clean.  I think that performing the PCR at a higher 
annealing temperature might help eliminate one of the contaminating bands.


Results
=======
.. datatable:: 20201016_optimize_ta_f12.xlsx

- I compared the results from the Sapphire (which better visualizes the minor 
  600 bp product) and the Geldoc (which doesn't over-saturate the 400 bp 
  product) to compare the different annealing temperatures.

- The best amplification is at :math:`T_a = \pu{56.6°C}`.

- The most specific amplification :math:`T_a = \pu{55.1°C}`, although for the 
  most part there is about 10% off-target activity for each annealing 
  temperature.

- Unexpectedly, the amount of off-target amplification increase at higher 
  annealing temperatures (although the difference is subtle).  I thought that 
  higher annealing temperatures would improve specificity.

Discussion
==========
- The best annealing temperature for this reaction is 55°C.

- Optimizing the annealing temperature did not significantly reduce the 
  off-target peak.


