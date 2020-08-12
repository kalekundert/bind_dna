*****************************************
Ligate linker-N with optimized conditions
*****************************************

Previously, I determined reasonable salt and temperature conditions for 
annealing, see :expt:`20191219_anneal_linker_n_with_salt` and 
:expt:`20200129_anneal_linker_n_with_temperature`.  These conditions gave rise 
to efficient ligation with the pseudo-linker, see
:expt:`20200129_ligate_with_5_phosphorylation`.  Here I want to test if the 
same protocol also gives efficient ligation for linker-N.

Methods
=======

Ligate Linker-N --- 2020/02/13
------------------------------
.. protocol:: 20200213_zap_anneal_ligate_formpage.txt

   See binder.

.. figure:: 20200213_ligate_o100.svg

.. datatable:: 20200213_ligate_o100.xlsx

   "Band 1" and "Band 2" refer to the uppermost and second uppermost bands on 
   the gel (both of which presumably represent ligated/annealed RNA/DNA).

- Ligation with linker-N is clearly much less efficient that it is with the 
  pseudo-linker.  I can only think of a few reasons why this might be:

   - The entropy of the flexible arm weakens binding of the linker to the mRNA.  
     If this is the case, the remedy would probably be to increase the amount 
     of complementarity, although I'm uneasy about making significant changes 
     to linker-N---the optimization could get really expensive.

     Note that this linker-N is missing the `CCTTG` arm meant to prime reverse 
     transcription, due to an error I made in ordering.  Perhaps this arm 
     really helps anchor the linker, even though it's not very long.  The 
     pseudo-linker does have this arm.

     Also note that the pseudo-linker has a 3' puromycin modification.  It may 
     have been smart to have ordered an internal Cy5 at the same position as 
     the 5-Me-dC phosphoramidite in linker-N, to mimic the mismatch there.  
     It's not a big deal, though.

   - My linker-N wasn't synthesized correctly.  I'm hesitant to blame Midland 
     CRC, because I've never had a company just send me the wrong thing, but 
     I've also never worked with this type of company before.  I'm also just a 
     little uneasy because I have no way to verify that linker-N is correct, 
     other than doing these reactions.  

     I'll be ordering more linker-N anyways, so if this is really the problem, 
     they'll have to make the same mistake again.  I can also ask them to use 
     MS for QC.  If I'm really worried about this, I could try ordering from 
     BEX instead (the company [Naimudden2016]_ used).  I could also try writing 
     to `Tai Kubo <mailto:tai.kubo@aist.go.jp>`_---the corresponding author on 
     the 2020 methods paper on cDNA display---and asking for some.

- I can't compare the upper and lower bands by gel densiometry, e.g. to compute 
  the ratio of linker-N that was ligated, because the lower bands are 
  saturated.  I can compare the upper bands, but I have to assume that I added 
  the same amount of linker-N to both lanes (which should be a pretty good 
  assumption).

  There is 2.3x more linker-N in the upper bands for the +ligase condition.  

- I don't know why there are two upper bands.  I expected only 1, which would 
  represent linker-N ligated (or possibly just annealed) to mRNA.  

  I wonder if one bands is ligated and the other is annealed but not ligated.  
  The idea would be that the denaturing conditions mostly prevent linker-N from 
  annealing, so ligated product is fully singled-stranded (despite having a lot 
  of complementarity).  This fully singled-stranded product might run 
  differently than annealed product.  If this were the case, though, I wouldn't 
  know which band is which, and I wouldn't know how to tell.  And in any case, 
  there clearly isn't that much ligated or annealed product.
  
- I don't know why linker-N seems to run slightly different in the ligation 
  reactions.  It's either a bit faster, or distributed a bit differently.  
  Unfortunately I kinda overran the gel; it would've been nicer to see 
  everything on the gel.

Ligate linker-N & click linkers --- 2020/08/10
----------------------------------------------
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
  
