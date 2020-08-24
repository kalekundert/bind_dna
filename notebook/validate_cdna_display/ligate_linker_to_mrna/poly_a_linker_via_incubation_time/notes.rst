*********************************
Poly-A linker via incubation time
*********************************

In this experiment, I want to test if I can improve ligation efficiency by 
incubation the reaction for different times at different temperatures.

Results
=======
2020/08/19:

.. figure:: 20200820_optimize_ligation_time_temp.svg

.. datatable:: 20200820_optimize_ligation_time_temp.xlsx

  The "Cy5 Only" columns were calculated from an image of the Cy5 channel taken 
  before GelGreen was added.  The "Cy5" columns were calculated from an image 
  of the same gel, after GelGreen staining.

- The reaction appears to be mostly complete after 10 min.

- The 25°C reactions went very slightly (~1%) further, but the difference 
  probably isn't significant.  I'll keep using 25°C because it's more 
  convenient.
  
- The RNA is noticeably degraded after the 24h incubation.  The non-uniform 
  banding of the degradation products suggests that this is due to RNase 
  contamination (as opposed to the RNA being hydrolyzed), and so it may be 
  possible to resolve this by adding RNase inhibitor.  But there doesn't seem 
  to be any benefit of a longer incubation time anyways.

- I inadvertently bleached the mRNA+linker and mRNA bands by scanning them 4-5 
  times in an effort to pick a good intensity threshold.  This can be seen by 
  comparing the yields calculated for the Cy5 channel before and after these 
  scans (see above, "Cy5 Only" vs "Cy5").  I was trying to get a more accurate 
  idea of saturation than I could get from the prescan, but it was a big 
  mistake to repeatedly scan only one of the two bands I wanted to compare.  

- Both the Cy5 and GelGreen channels indicate about 55-60% yield, which is what 
  I've seen in previous experiments.  There was a bigger disagreement for the 
  "0m" reaction, but this may in part be because the signal in the GelGreen 
  channel was pretty dim for this sample.

- The "0m" reaction didn't completely eliminate ligation.  That's probably 
  because this was more like a 1-2m reaction, since it took me some time to: 
  (i) add the annealed mRNA/linker to the ligation master mix, (ii) mix it, 
  (iii) take an aliquot for the 0m reaction, (iv) add EDTA, (v) mix it, (vi) 
  out it in the freezer, and (vii) freeze.  If I were to do this again, I'd 
  probably flash freeze the sample, since I really don't know how long that 
  took.  I might also not use the master mix for the 0m reaction (although then 
  I have to confront the fact that I can't accurately pipet very small volumes 
  of ligase).

The Cy5 and GelGreen channels agree pretty well with each 
  
Because the linker band was not included in the region being scanned, and was 
therefore not bleached, comparisons between the was a particular This wasn't a 
big problem for interpreting the GelGreen channel, because 

- I bleached Cy5 by trying to optimize the intensity setting on the laser 
  scanner.  This is evident by comparing the "Cy5" and "Cy5 Only" yields.

- The RNA is significantly degraded after a 24h incubation at either 16°C or 
  25°C.  It may be possible to prevent this degradation by adding RNase 
  inhibitor.

