*********************
Block BAL-31 via TelN
*********************

It would be convenient to be able to block BAL-31 activity, for reasons 
discussed in :expt:`20190419_validate_exonuclease_assay`.  The TelN enzyme is a 
promising way to do this because it makes cuts DNA in a way that leaves 
covalently capped ends.  Since BAL-31 can only cut at DNA ends, and TelN leaves 
no ends, it is reasonably to think that TelN might block BAL-31 activity.

Results
=======

TelN & BAL-31 --- 2019/11/15
----------------------------
.. protocol:: 20191115_phenol_chloroform_extraction.txt 20191115_pcr.txt 

   See binder for TelN digestion (2019/11/12), BAL-31 digestion, E-gel 
   parameters, and discussion of how to clean up the reaction.

.. figure:: 20191115_bal31_blunt_vs_capped.svg

   Note that I enhanced the contrast of the lanes with DNA.

Discussion
==========
- There is no evidence that TalN blocks BAL-31 activity.  This is surprising to 
  me, so I might want to repeat this experiment, but the results are pretty 
  convincing.
