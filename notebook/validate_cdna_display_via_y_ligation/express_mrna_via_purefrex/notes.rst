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

PUREfrex 1.0, mWasabi
=====================

PUREfrex 1.0, FLAG
==================

Discussion
==========
- The optimal mWasabi mRNA concentration for a PUREfrex reaction is ≈500 nM.  1 
  µM mRNA gives only slightly higher expression.  Adding RNase inhibitor does 
  not affect this concentration.

- For some reason, the PUREfrex kit is the only one that does not give a 
  homogeneous product for mWasabi.
