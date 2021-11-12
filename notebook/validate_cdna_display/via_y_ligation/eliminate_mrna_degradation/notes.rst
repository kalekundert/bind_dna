**************************
Eliminate mRNA degradation
**************************

In a number of experiments, including :expt:`14`, I have seen the mRNA degrade 
significantly during the ligation reaction.  I suspext that this is due to 
RNase contamination.  My goal for this experiment is to:

- See if RNase inhibitor helps.

- Identify the reagent causing the degradation.

- If possible, acquire that reagent without contamination.

Results
=======

Reagent dropout --- 2020/01/27
------------------------------
I setup annealing and ligation reactions, leaving one reagent out of each 
reaction.  I also setup a control with mRNA (f11) and linker (o93) mixed in 
water just before running the gel, so check how well the loading dye is 
denaturing secondary structure.

.. protocol:: 20200127_reagent_dropout_page.txt

   See binder.  I used 500 mM NaCl instead of 1x PBS by mistake.  A few other 
   details in the printed protocol were either inaccurate or missing.

.. figure:: 20200127_ligation_reagent_dropout.svg

- None of the lanes exhibit significant degradation.

  It may be that the RNA degradation I saw it the past was just to due one-off 
  contamination of those reactions, not my stocks.  Or it could have been due 
  to the PBS, which I didn't use this time (by mistake).  The pH of PBS is 7.4, 
  which should be ok with RNA.

  In any case, it occurred to me that the better way to test for RNase 
  contamination would be to mix each reagent individually with the mRNA, not to 
  do a dropout like this.  If there were multiple contaminated reagents, the 
  dropout approach wouldn't be informative.

- The urea loading buffer I've been using doesn't seem to be completely 
  breaking the base pairing between the mRNA and the linker.  Note that the 
  leftmost lane has a bright yellow band, despite not being annealed or 
  ligated.  Similarly, the −ligase reaction also has a bright yellow band.  
  Fitzy suggested that a formamide loading buffer might denature the RNA more 
  completely.

  .. note::

     The NEB RNA loading dye (B0363), which I think came with my ladder, is 
     formamide-based.  So I can try this right away.

  That said, the brightest linker bands are mostly in the lanes that either 
  don't have dialysis or that you'd expect to have to poorest ligation (the 
  exception being the −ligase lane, which has only a faint linker band).  This 
  suggests that annealed linker can survive dialysis, then be partially removed 
  by the loading buffer.

  This data also hints that the intense annealing protocol I'm doing may not 
  really be necessary.  The mRNA and the linker seem to anneal pretty easily, 
  any many annealing protocols in the literature are not as intense.  I can 
  probably make mine simpler/shorter, but I'd rather optimize that once the 
  ligation is pretty much working.

- Next time, don't do dialysis.  Diluting is faster, and would make the linker 
  band more interpretable.
