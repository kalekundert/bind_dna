******************************
Linker-N, optimized conditions
******************************

Previously, I determined reasonable salt and temperature conditions for 
annealing, see :expt:`12` and :expt:`13`.  These conditions gave rise to 
efficient ligation with the pseudo-linker, see :expt:`16`.  Here I want to test 
if the same protocol also gives efficient ligation for linker-N.

Methods
=======
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

