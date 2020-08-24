**************************************
Poly-A linker via optimized conditions
**************************************

Previously, I determined reasonable salt and temperature conditions for 
annealing, see :expt:`12` and :expt:`13`.  These conditions gave rise to 
efficient ligation with the pseudo-linker, see :expt:`16`.  Here I want to test 
if the same protocol also gives efficient ligation for linker-N.

Results
=======

Formamide buffer --- 2020/08/08
-------------------------------
.. protocol:: 20200808_rnasezap_anneal_mrna_linker_ligate_linker_n_gel.txt

.. figure:: 20200810_ligate_linkers_post_stain.svg

- I can't say for sure that the ligation reaction worked, because the formamide 
  denaturation wasn't effective.  I should look to see how [Naimudden2016] 
  demonstrated that their ligation was successful, if they did.
  
  [Wang2014]_ describes various approaches for denaturing dsDNA into ssDNA, and 
  it seems like this is actually not a trivial thing to do (if you want the DNA 
  to not immediately renature).  For the purposes of this experiment, the best 
  approach seems to be adding DMSO.  DMSO and formamide are believed to work in 
  similar ways, although DMSO is a slightly stronger denaturant.  Both solvents 
  are more effective when used at higher concentrations, suggesting that I may 
  be able to get more mileage out of formamide by just using more of it (i.e.  
  for this experiment I used a 1:1 formamide:sample ratio; I could try 
  increasing that to 4:1 or 9:1).

  I can't find any examples of anyone using DMSO in a loading buffer, let alone 
  using a loading buffer composed of >95% DMSO.  So I won't be surprised if 
  this doesn't work well for some reason.  DMSO (1.10 g/mL) and formamide (1.13 
  g/mL) have similar densities, so at least the sample should sink in the well 
  with either solvent.
  
  Alternatively, I could try incubating at 95°C after adding a significant 
  excess of free oligo with the Y-tag sequence (o194).  DNA does denature at 
  high temperatures, but renatures rapidly once the temperature is cooled.  
  This approach would give the DNA a way to renature that would still allow me 
  to tell whether or not there is a covalent attachment between the linker and 
  the mRNA.  One drawback is that the free oligo would appear on the gel, and 
  would be very bright (due to being present in significant excess).  It should 
  run faster than the free linker (i.e. it should run just like o93), so I 
  could try running it off the bottom, but that would be delicate and would 
  probably take me a few tries to get right.

- The amount of annealed DNA seems significantly worse in this experiment than 
  it did in :expt:`51`.  I wonder if this is due to differing reaction 
  conditions.

  In the aforementioned experiment, I didn't use any salts: just water and the 
  oligos.  In this experiment, the annealing step has 137 mM NaCl and the 
  ligation step has 10 mM MgCl₂.  The salt in the annealing step may not be 
  relevant, because its gone by the time the samples are being prepared for 
  electrophoresis, but the MgCl₂ might be.  Millimolar concentrations of 
  divalent salts are known to `stabilize oligo duplexes 
  <https://www.idtdna.com/pages/education/decoded/article/understanding-melting-temperature-(t-sub-m-sub-)>`__.

- I might try doing a similar experiment with excess linker, so I can see if I 
  can drive the reaction to completion.  I'd rather not use a huge excess of 
  linker, because it is expensive, but I think it's worth experimenting with.

- The linkers with the puromycin arm (o100, o130, o129, o128) visibly shift the 
  mRNA, while the linkers without it do not (o93, o127) don't.  This makes 
  sense, and it nice to see.

- linker-N appears to have been ruined as a result of being left at room 
  temperature for an indeterminate period of time during the shutdown.  There 
  are a couple signs that things were wrong:

  - When I went to resuspend the oligo, it was already in a small volume of 
    liquid.  Presumably it was shipped lyophilized, but maybe it is hygroscopic 
    enough to pull water out of the atmosphere.  In any case, being dissolved 
    at room temperature probably left in even more prone to degradation.

  - The oligo is a darker shade of blue than my other oligos labeled with Cy5.  
    In combination with the fact that it doesn't appear in the gel, I think 
    it's fair to say that its Cy5 has been inactivated.

- The clicked linkers (o129, o128) anneal/ligate noticeably better than 
  linker-N − RT (o100).  Presumably this is because the RT-arm provides 
  important complementarity to the Y-tag sequence.  It would've been nice to 
  compare to linker-N with the RT-arm (o130), but that oligo seems to have been 
  destroyed as discussed above.

- It's interesting that the click reaction seems to proceed a little further 
  (~95% complete rather than ~90% complete) after the ligation reaction.  

  .. datatable:: percent_unclicked_after_ligation.xlsx
  
  Assuming this isn't some artifact, the difference could be due to the various 
  incubation steps at elevated temperatures, specifically:
  
  - 65°C for 10 min (to denature the ligase)
  - 70°C for 3 min (prior to running the gel)
    
  I looked briefly to see if there were reports of Cu-free click working better 
  at elevated temperatures, and while I didn't find anything really conclusive, 
  I did find one person recommending 37°C.  I think it would be reasonable to 
  experiment with different incubation temperatures (perhaps the next time I 
  buy more oligo; I don't know if I have enough at the moment).
  
Competing oligo --- 2020/08/10
------------------------------
See :expt:`50` for more details about this experiment.

.. figure:: 20200814_titrate_o194_saturated.svg
