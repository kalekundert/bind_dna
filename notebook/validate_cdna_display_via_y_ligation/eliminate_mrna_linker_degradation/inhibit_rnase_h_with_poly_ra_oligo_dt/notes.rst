*************************************
Inhibit RNase H with poly-rA/oligo-dT
*************************************

[Walder1988]_ shows that the RNase H present in cell lysates is responsible for 
degrading mRNA targeted by anti-sense oligos.  As part of this work, the 
authors use poly-rA/oligo-dT to inhibit RNase H activity.  I want to see if 
this approach works for me, although I'm worried about the potential for 
poly-rA/oligo-dT to interfere with reverse transcription.

Considerations
==============

How much poly-rA/oligo-dT?
--------------------------
- [Walder1988]_ used a final concentration of 8 A260 U/mL, but also saw pretty 
  good inhibition at 4 U/mL.

- I assume that the ratio of decoy substrate to mRNA is what's really relevant, 
  so I need to know how much mRNA was used by [Walder1988]_.

  - [Walder1988]_ used total mRNA extracted from mouse reticulocyte cells, at a 
    concentration of 20 µg/mL.

  - >95% of this mRNA is α/β globin [Walder1988]_, [Kabat1977]_.

  - α/β globin mRNA is 9S [Kabat1977]_.

  - Assuming S units are comparable between mRNA and rRNA, 9S would be about 1 
    kb.  This at least seems reasonable.

    https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-isolation/general-articles/ribosomal-rna-sizes.html

  - ssRNA MW ≈ 320 × N = 320,000 Da

  - This corresponds to a rough concentration of 60 nM.  This is very 
    approximate, but I wasn't able to find any better ways to estimate the 
    average length of the mRNA is the [Walder1988]_ experiments.

- Based on the results for :expt:`18`, I add mRNA to a final concentration of 
  160 nM in my PURExpress reactions.

- I consider these concentrations comparable, especially given how 
  approximate the [Walder1988]_ concentration is.

- PURExpress probably has less RNase H activity that rabbit reticulocyte 
  lysate, in which case I would need less inhibitor.

- All things considered, I think it makes sense to start with 8 A260 U/mL, and 
  if that works, to try some dilutions.

- I ordered 5 A260 U, so I'll probably want a stock of 80 U/mL, which would be 
  62.5 µL.  I'd use ≈1 µL per reaction.

Results
=======

2021/02/15
----------
.. protocol:: 20210215_protect_mrna_with_ra_dt.txt

.. figure:: 20210215_protect_mrna_with_ra_dt.svg

.. datatable:: 20210215_protect_mrna_with_ra_dt.xlsx

Observations:

- The poly-rA/oligo-dT seems to have a minor positive effect on mRNA stability, 
  and a minor negative effect on protein expression.  These differences are 
  hard to see in the gel itself, but are apparent in the densiometry analysis.

- I see once again that the "cleaved" Cy5 band runs a bit slower than the o129 
  band.  I think this is consistent with RNase H contamination.  See :expt:`65` 
  for a discussion of this.

- The GFP signal was very faint, probably for the same reasons discussed in 
  :expt:`62` (i.e. short reaction time).

- About 90% of the linker in this f89 prep is ligated to mRNA (data not shown).

Conclusions:

- It would've been nice to include an f85 control.  I'm curious if f85 would 
  give a different level of mWasabi expression than f89.

- I'm still not sure what going on.  On one hand, if there was RNase H 
  contamination, I would've expected poly-rA/oligo-dT to have a bigger 
  inhibitory effect.  On the other hand, if there wasn't RNase H contamination, 
  I wouldn't have expected poly-rA/oligo-dT to have any effect at all.  The 
  weak effect that I observe could go either way.

- I wonder if it would be informative to do similar experiments with decoy 
  ssDNA, dsDNA, ssRNA, and dsRNA.
