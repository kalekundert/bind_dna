********
Via salt
********
As discussed in :expt:`53`, salt is widely used in oligo annealing reactions to 
help shield the charge of the phosphate backbone.  However, [Naimudden2016]_ 
doesn't clearly specify the buffer used to anneal the oligos.  Therefore, this 
is a parameter that I need to optimize myself.

Results
=======

Removing salt --- 2019/12/19
----------------------------
A significant issue with the annealing reaction is that salt is inhibitory to 
the phosphorylation step, and may also be inhibitory to the ligation step.  See 
more discussion in :expt:`19`.  To find the best way to address this issue, I'm 
trying 4 different strategies for removing/diluting the salt:

- A: Pre-phosphorylate the linker; add 50-150 mM NaCl to the annealing 
  reaction; don't bother removing it before the ligation reaction.

- B: Anneal in the smallest possible volume of PBS; drop dialysis to remove 
  salt; ligate in T4 DNA ligase buffer.

- C: Anneal in the smallest possible volume of PBS; dilute 10x for ligation 
  reaction.

- D: Anneal and ligate in T4 ligase buffer (i.e. without salt).  I think this 
  is what [Naimudden2016]_ did.

Note that I'm using 5 pmol of mRNA and linker-N in each reaction instead of the 
50 pmol called for by [Naimudden2016]_.  I'm scaling down because I don't have 
enough RNA to do 4 reactions with 50 pmol.  And anyways, since [Naimudden2016]_ 
didn't specify a reaction volume, I can just scale up whatever volume works for 
me when I go back to doing 50 pmol reactions.

.. figure:: 20191219_optimize_linker_n_ligation.svg

- I really don't know why the linker-N control seems to have both mRNA and 
  ligated mRNA.  It's concerning; I'm pretty certain I didn't mix up samples (I 
  prepared the controls just at the end, and while the other reactions were in 
  the thermocycler) and I really don't want my linker-N to be contaminated.

  .. update:: 2020/01/29

      I haven't seen this mixture in subsequent experiments, so just messed up 
      some pipetting; I didn't contaminate the stock.  I should think about 
      making aliquots, though.

- All of the lanes should have the same amount of linker-N, except the lanes 
  that were diluted (C1 and C2 intentionally, B1 by mistake).  But the A and D 
  lanes seem to have very little mRNA.

  Actually, this is an illusion.  When I subtract the background and add up all 
  the pixels in each lane, I see that A == D >= B2 > B1 > C1 == C2.  This makes 
  sense: B1 should have 0.65x (check that number) the amount of RNA in B2, due 
  to the pipetting error I made, and this is about what I see.  C1 and C2 
  should have 0.25x the amount of RNA in the other lanes, because they were 40 
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

NaCl, PBS, [Co(NH₃)₆]Cl₃ --- 2020/01/29
---------------------------------------
[Kitamura2002]_ uses [Co(NH₃)₆]Cl₃ in the ligation buffer.  According to 
`Wikipedia 
<https://en.wikipedia.org/wiki/Hexamminecobalt(III)_chloride#Uses>`_, 
"[Co(NH₃)₆]³⁺ is an unusual example of a water-soluble trivalent metal complex 
and is of utility for charge-shielding applications such as the stabilization 
of highly negatively charged complexes, such as interactions with and between 
nucleic acids."  So there is reason to think that Co³⁺ may be more effective 
for annealing than Na⁺.

- For this experiment, I'll use the "95°C→25°C over 1h" thermocycler protocol.  
  I haven't optimized the temperature gradient yet, but this one seems to be 
  very in line with other protocols and I don't expect that it will cause 
  problems.

- I decided to use linker-N rather than the pseudo-linker (o93) because it 
  seems from :expt:`11` that linker-N may not anneal as tightly to the mRNA as 
  the pseudo-linker does.  If there's a significant difference between the 
  oligos, I care more about the behavior of linker-N, and I'd rather see the 
  effect of salts on the weaker binding oligo.

.. protocol::

   See binder, 2020/01/29

.. figure:: 20200129_anneal_na_pbs_co.svg

.. datatable:: 20200129_anneal_na_pbs_co.xlsx

- I expected a more significant fraction of the linker to be annealed.

   - I wonder somewhat if I'm adding too much linker to the reaction.  I'd 
     really have to be off by a factor of 100 or something for that to fully 
     explain the relative faintness of the upper bands.  And the error would 
     probably have to be in the mRNA, since both linker-N and o93 seem 
     similarly concentrated (and I probably didn't make the same 100x dilution 
     error twice).

   - The relatively poor annealing could also just be because o100 is missing 
     its RT-arm.  The poor annealing could then possibly explain the poor 
     ligation efficiency.  I definitely need to order the right oligo.

- I think the high MW bands (~800) are mRNA dimers.

   - High salt is consistent with more base pairing in nucleic acids.

   - Maybe I can ask Vienna what it thinks the dimer would be.

   - I don't want mRNA to be ligated together, but so far I haven't seen that 
     in any of my attempted ligation reactions.

   - In the gel densiometry results above, I combined green pixels from both 
     bands, since any green outside the lowest bands must represent linker-N 
     annealed to mRNA.

- PBS seems to work well.

   - It's interesting that PBS looks quite different than 137 mM NaCl.  Perhaps 
     the relatively low concentration of Mg²⁺ has an outsized effect.

- 500 mM NaCl, despite not have the most annealed pixels, is the only condition 
  that has a discernible mRNA band corresponding to its linker-N band (not 
  counting the ≈800nt bands).  It might be worth doing the whole ligation 
  reaction with 500 mM NaCl, to see how well it works.

- It's interesting that the reaction without salt *and* without annealing seems 
  to do pretty well.  Especially since the reaction without salt *but* with 
  annealing performs much worse.  Is there something unexpected going on?  Or 
  is this assay maybe just noisy?

   - Also note that both the "f11 only" and "no salt, no temperature" controls 
     have faint 800 nt bands in the GelRed channel.  This really makes it seem 
     like theres something about no salt and no temperature that allows for 
     annealing (since it's more than just one reaction).

   - Maybe I should try with salt and without temperature.

- Cobalt seems to destroy the mRNA.  This is the second time I've seen this, so 
  I'm definitely starting to think that Co catalyzes the cleavage of the RNA 
  backbone or something.

Discussion
==========
- I'm tentatively planning to use PBS for future experiments.
