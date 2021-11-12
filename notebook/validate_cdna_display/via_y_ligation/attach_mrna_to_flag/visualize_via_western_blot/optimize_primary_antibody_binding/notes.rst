*********************************
Optimize primary antibody binding
*********************************

Before doing any experiments on unknown samples, I want to identify good 
parameters for visualizing known samples.  Specifically, I want to optimize the 
concentration and incubation time for the primary antibody.

.. protocol:: 20210923_optimize_western.pdf 20210923_optimize_western.txt

.. figure:: 20210924_optimize_western.svg

- The ladder is visible, but the FLAG peptide is not.  This rules out problems 
  with the imager and (for the most part) the transfer.  Remaining hypotheses 
  of why I don't see FLAG:

  - I didn't load enough FLAG.
  - Either the primary or the secondary didn't work.
  - The FLAG blew through the membrane.
  - Too much washing.

- I forgot to mark which side of the membrane was facing the gel, so I tried 
  imaging both sides.  I didn't see any significant difference between the two.



