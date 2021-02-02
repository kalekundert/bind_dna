**************************
Test extracts via Wang2017
**************************

Now that I've shown that I can detect RNase H activity using the [Wang2017]_ 
assay, I want to see if I can detect RNase H activity in any of my *in vitro* 
translation systems.

Results
=======

2020/12/16
----------
.. protocol:: 20201216_detect_rnase_h.txt

.. figure:: 20201216_test_extracts_for_rnase_h_fits.svg
.. figure:: 20201216_test_extracts_for_rnase_h_slopes.svg

To do:
- Should repeat with the following controls:
  - −DNAzyme, to really see what the baseline fluorescence looks like.
  - −beacon, −DNAzyme same idea.
  - o212, to see if signal is detectable.
  - o211, really just checking that the assay works.

2021/01/13
----------
.. protocol:: 20210113_detect_rnase_h.txt

.. figure:: 20210113_test_extracts_for_rnase_h_fits.svg
.. figure:: 20210113_test_extracts_for_rnase_h_slopes.svg

Observations:

- The lysates have very high fluorescent backgrounds, and the linear fits for 
  these curves do not match the data well at all.

- NEBExpress has low fluorescent background, but a fairly significant slope in 
  all conditions.  This may indicate that something in NEBExpress can cleave 
  the molecular beacon (o213) on it's own.

- The slopes in the water controls are greater than they have been previously 
  (~4.1 and 2.3 this time; 2.1 and 1.3 on 2020/12/16; 1.4 and 0.7 on 
  2020/12/15).  It's also interesting that the o210 (RNA-blocked DNAzyme) 
  control consistently has a slope ~2x greater than the o211 control 
  (DNA-blocked DNAzyme).  This might reflect the intrinsic difference in 
  stability between RNA and DNA.

- In the PURExpress conditions, the positive control (o212) has the same signal 
  as the negative control (o211), even as the experimental sample (o210) has 
  much higher signal than either.  I can't make sense of this, especially 
  considering that the positive and negative controls worked well in the water 
  conditions.

Discussion
==========
- The results for the lysates are meaningless and should be ignored.

- I won't be able to use assays with fluorescent readouts to measure RNase H 
  activity in lysates.  [Wang2007]_ claimed to measure RNase H activity in 
  lysates, and got around the background fluorescence problem diluting the 
  lysate 100-fold.  I don't think the assay would be sensitive enough to detect 
  activity if I did the same, though.

  I might be able to measure RNase H activity in these samples using a qPCR 
  based assay.

- There may be RNase H activity in PURExpress, but I'm hesitant to make any 
  strong conclusions because the positive control doesn't make sense.  That 
  said, in both PURExpress replicates the o210 well had a significantly greater 
  slope that the o211 sample.
