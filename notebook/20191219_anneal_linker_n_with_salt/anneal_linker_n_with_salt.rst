*************************
Anneal linker-N with salt
*************************
Trying 4 different strategies for linking to see which gives the best 
efficiency:

- A: Pre-phosphorylate the linker, add 50-150 mM NaCl to the annealing 
  reaction, don't bother removing it before the ligation reaction.

- B: Anneal in the smallest possible volume of PBS, drop dialysis to remove 
  salt, ligate in T4 DNA ligase buffer.

- C: Anneal in the smallest possible volume of PBS, dilute 10x for ligation 
  reaction.

- D: Anneal and ligate in T4 ligase buffer (i.e. without salt).  I think this 
  is what [Naimudden2016]_ did.

Results
=======
Note that I'm using 5 pmol of mRNA and linker-N in each reaction instead of the 
50 pmol called for by [Naimudden2016]_.  I'm scaling down because I don't have 
enough RNA to do 4 reactions with 50 pmol.  And anyways, since [Naimudden2016]_ 
didn't specify a reaction volume, I can just scale up whatever volume works for 
me when I go back to doing 50 pmol reactions.

.. figure:: xxx

- I really don't know why the linker-N control seems to have both mRNA and 
  ligated mRNA.  It's concerning; I'm pretty certain I didn't mix up samples (I 
  prepared the controls just at the end, and while the other reactions were in 
  the thermocycler) and I really don't want my linker-N to be contaminated.

- All of the lanes should have the same amount of linker-N, except the lanes 
  that were diluted (C1 and C2 intentionally, B1 by mistake).  But the A and D 
  lanes seem to have very little mRNA.

  Actually, this is an illusion.  When I subtract the background and add up all 
  the pixels in each lane, I see that A == D >= B2 > B1 > C1 == C2.  This makes 
  sense: B1 should have 0.65x (check that number) the amount of RNA in B2, due 
  to the pipetting error I made, and this is about what I see.  C1 and C2 
  should have 0.25x the amount of DNA in the other lanes, because they were 40 
  µL reactions rather than 10 µL, and this is again what I see.  Note that I 
  can't compare to the mRNA control lane, because it's saturated.

  .. datatable:: 20191219_lane_densiometry.csv

      Values are for the entire lane (except the linker-N bands), i.e. all mRNA 
      regardless of how diffuse it is.

   .. note::

      Add corrected columns to this table.

  So all the lanes have about the same amount of RNA, but the RNA in lanes A 
  and D is much more diffuse for some reason.

- I don't know why all of the ligation products (especially A and D) are so 
  much more diffuse than the initial mRNA (which is quite clean).  I could do 
  no-annealing/no-ligation controls to find out if it's related to either of 
  those steps.  Maybe I need to dilute the reaction to prevent intermolecular 
  ligation.

- The B and C lanes (once the C lanes are multiplied by 4 to account for their 
  dilution) gave the most ligated product.  Especially since D gave no product, 
  this strongly suggests that salt in the annealing step is essential.  

- The B lanes have much less free linker-N than the other lanes, even the 10x 
  diluted lanes.  Part of this is because more of the linker was ligated to the 
  mRNA, but a bigger effect (I think) is that un-annealed linker-N is being 
  lost during dialysis.  Note that the B lanes have even less un-annealed 
  linker-N than the C lanes, which are diluted 4x.

- Maybe run the gel for 70+ minutes next time.  This would run un-ligated 
  linker off the bottom, but would give better separation between the mRNA 
  bands.

Optimize dialysis --- Planning
------------------------------
Does dialysis really help?  It seemed to be the best method in the 12/19 
experiment, but I want to run the following controls to better understand why: 

- No salt: Is dialysis helpful on its own?

- No dialysis: Is it necessary to get rid of salt?

- No ligase: 

- 
