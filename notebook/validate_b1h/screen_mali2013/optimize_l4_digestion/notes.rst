*********************
Optimize l4 digestion
*********************

2022/09/22:

When gel-purifying f191 for the first time, I noticed 3 bands around the 
expected MW.  I suspect that these bands represent partial digestion products.  
To confirm that, I want to vary the amount of time that the digestion reaction 
proceeds for.  (The protocol calls for 5-15 min, and I initially did 15 min.)

Considerations
==============
How many timepoints?

- The `NEB restriction digest troubleshooting guide`__ mentions that, if 
  necessary to increase incubation times, 1-2h is typically sufficient.

__ https://www.neb.com/tools-and-resources/troubleshooting-guides/restriction-enzyme-troubleshooting-guide

- Want to see a progression, so want a couple different times:

  - 0m
  - 5m
  - 15m
  - 30m
  - 1h
  - 2h

How much material?

- The three bands are faint in my current gel, which loads the entire PCR 
  reaction in a single lane.  

- I have a PCR reaction that I ran for 35 cycles.  I could use that, but will 
  the extra cycles affect anything?  Once the PCR reaction runs out of primers, 
  you can get library members annealing to each other, leading to dsDNA with 
  bulges.  These bulges would be adjacent to the BbsI site, so perhaps they 
  would matter.

- The alternative is to run 6 25 µL reactions.  I can do that; I have tons of 
  reagents.

Results
=======
.. protocol:: 20220923_optimize_l4_digestion.pdf 20220923_optimize_l4_digestion.txt

.. figure:: 20220923_bbsi_timecourse.svg

.. datatable:: 20220923_bbsi_timecourse.xlsx

- The reaction is complete after 5 min.

- Only ≈90% of the DNA is fully digested.

  - The band at ≈260 bp is probably due to partial digestion.  These products 
    would be 264 bp, which is pretty close to what we see.

  - This probably isn't an issue with the enzyme, since it's clearly active and 
    finishes the reaction in less than 5 min.

  - There could be mutations in the restriction sites that prevent cleavage.  
    That would imply a pretty high error rate, though.

    - l4 is 264 bp, of which 12 bp form the two restriction sites.
    - A 10% chance that one of those sites has a mutation (given that 90% of 
      the DNA is fully digested) would roughly be a 1% chance that each 
      position has an error, or a 1/100 error rate.
    - I couldn't find a published error rate for oPools, but 1/100 is way too 
      high.

- 90% yield should be good enough, though.
