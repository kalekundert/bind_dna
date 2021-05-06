*************************
Express mRNA via PUREfrex
*************************

Before I an do experiments with PUREfrex, I need to determine the optimal 
amount of mRNA to use in reactions.  My plan is to closely follow the protocol 
from :expt:`18` (2020/02/26).

PUREfrex 2.0, mWasabi --- 2021/02/19
====================================
.. protocol:: 20210219_optimize_purefrex2_s30.txt

.. figure:: 20210219_optimize_purefrex2_s30.svg

.. datatable:: 20210219_optimize_purefrex2.xlsx

Observations:

- There are two green fluorescent bands.

  By comparison to the Promega S30 gel, the upper band is probably full-length 
  mWasabi.  Note however that the S30 band runs slower than either PUREfrex 
  band, so it's possible that neither band is full-length.

  I can't really say exactly how big the mWasabi fragments are, because mWasabi 
  always runs slower than the comparable ladder bands.  I initially thought 
  that this might be due to the fluorophore rearrangement creating an 
  unbreakable "crosslink", but the reaction just cyclizes two adjacent 
  residues.  So I don't know why there's a discrepancy between the ladder and 
  mWasabi.

  Interestingly, the mWasabi expressed using PURExpress seems to run slightly 
  differently than any of the bands in the above gels.  (I would put it between 
  the two PUREfrex bands, roughly.)  The NEBexpress bands seem to run about the 
  same as the S30 bands (which makes sense, given that both are S30 extracts).  
  So maybe something about the IVTT mix affects how mWasabi runs.

  If I had to guess, I'd say that the lower band is the product of an internal 
  start site.  But if that were the case, I'd expect the upper band to become 
  the major product at lower mRNA concentrations.  Instead, the products appear 
  to be in about the same ratio at each mRNA concentration.

PUREfrex 2.0 + RNase inhibitor, mWasabi --- 2021/02/23
======================================================
.. protocol:: 20210223_optimize_purefrex2_rnase.txt

.. figure:: 20210223_optimize_purefrex2_rnase.svg

Observations:

- Ratio of small to large mWasabi bands is about equal in all conditions.

- Relative expression levels are comparable to the −RNase inhibitor reactions 
  from 2021/02/19.

PUREfrex 1.0, mWasabi --- 2021/04/14
====================================
.. protocol:: 20200414_optimize_purefrex1_gfp.txt

.. figure:: 20210414_optimize_purefrex1_gfp.svg

Observations:

- Unlike with PUREfrex 2.0, with PUREfrex 1.0 there is only a single GFP band.

- Looks like much lower GFP expression relative to PUREfrex 2.0, based on the 
  fact intensity of the mWasabi bands relative to the ladder.  I'd need to 
  compare both in the same experiment to know for sure, though.  As noted by 
  [Reyes2021]_, high expression is not necessarily desirable for DNA display.

Conclusions:

- 1000 nM is the optimal mRNA concentration for mWasabi with PUREfrex 1.0.

PUREfrex 1.0, FLAG --- 2021/04/07
=================================
.. protocol:: 20210407_optimize_purefrex1_flag.txt

.. figure:: 20210407_optimize_purefrex1_flag.svg

Observations:

- I'm not sure whether or not I can see the FLAG peptide product.  There's a 
  band at 3 kDa that seems to get brighter with increased template 
  concentration, but the band is present even in the −template control.  Also, 
  the FLAG peptide is 2.2 kDa, so the band isn't quite in the right spot.  
  Maybe peptides don't run exactly at their MW, though.

- The loading dye is very clearly visible, which annoys me.  I'm going to try 
  avoid this problem by making a crystal violet loading dye next time 
  [Tice1991]_.

Conclusions:

- I need to try adding RNase to see if that helps the visualization.  I won't 
  be able to do that when I also want to visualize mRNA, but it could still be 
  useful for experiments like this (or once I have cDNA).

PUREfrex 1.0, Zif286
====================
.. protocol:: 20210505_optimize_purefrex1_zif.txt

.. figure:: 20210505_optimize_purefrex1_zif.svg

Observations:

- The optimal mRNA concentration is about 1000 nM.  This is in line with all of 
  my other experiments.

- The Zif268 is discernible with LysGreen and tricine PAGE.

  - This is a big improvement over Coomassie staining with bis/tris/MES PAGE 
    (:expt:`18`), which doesn't detect Zif268 at all.

  - The bands are still very faint (you can see how high the contrast is turned 
    up).  Part of the problem is just that PUREfrex 1.0 has relatively low 
    expression.  I'd probably see a stronger band with PUREfrex 2.0.

Discussion
==========
- For every condition, 1 µM seems to be the optimal mRNA concentration.  
  However, 500 nM often gives comparable expression, and so might be preferred 
  in the interest of conserving material.

- Adding RNase inhibitor does not affect the optimal mRNA concentration.

- PUREfrex 2.0 does not give a homogeneous product for mWasabi.  Every other 
  kit I've tried—including PUREfrex 1.0—does.
