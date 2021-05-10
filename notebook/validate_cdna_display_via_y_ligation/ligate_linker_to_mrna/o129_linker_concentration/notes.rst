**************************
o129, linker concentration
**************************
Although the annealing and ligation reactions seem to work with the poly-A 
linker (o129), only about half of the mRNA is modified.  An obvious way to try 
modifying more of the mRNA is to use more linker.

Considerations
==============
[Naimudden2016]_ reports that other groups have used between 2–200x excess of 
linker for Y-ligation (refs 13, 20, 26).  In my own survey of :expt:`T4 RNA 
ligase protocols <-1>`, I found that most protocols recommend no excess of 
linker (i.e. 1:1 RNA:linker), but one (from NEB) recommends a 2-10x excess.

Disadvantages of using excess linker:

- Needs to be removed, because any extra puromycin will interfere with the 
  protein expression reaction.  However, the 100 kDa MWCO spin column should do 
  this pretty effectively (o129 is about 20 kDa).

- The linker is one of the most expensive reagents.

Considering the  above, I want to limit myself to a fairly narrow range of 
linker concentrations.  Specifically, I'll start by testing:

- 1x
- 2x
- 4x
- 8x

Results
=======
2020/09/03:

.. protocol:: 20200903_optimize_linker_conc.txt

.. figure:: 20200903_optimize_linker_conc.svg

.. datatable:: 20200903_optimize_linker_conc.xlsx

- Adding excess linker has does not increase the yield of the ligation 
  reaction.  In fact, if anything it seems to diminish it.

- My serial dilution wasn't perfectly accurate.  It looks like I diluted by 
  closer to 2.5x on each dilution.  I recalculated the amount of excess linker 
  in each reaction assuming that the most concentrated linker (80 µM) was the 
  most accurate.

- It's interesting that my yield didn't even decrease when I only added 0.5x 
  linker.  This indicates that it's something about roughly half of the mRNA 
  that is preventing the reaction from going to completion.  Possible 
  explanations:

  - T7 polymerase is adding extra nucleotides to the 3' end of the mRNA 
    [Gholamalipour2018]_.  I can't find the reference anymore, but I remember 
    reading that the 5' and 3' sequence is important for Y-ligation.  If this 
    is the case, it's reasonable to think that these extra nucleotides could 
    inhibit ligation.  Note that the high-yield conditions of the T7 HiScribe 
    reaction increase the probability of there being extra nucleotides 
    [Gholamalipour2018]_.  Possible solutions:

    - Use a cis-ribozyme.  See :expt:`64`.

    - Cleave with oligo-directed RNase H or tRNA-directed RNase P.  
      Cis-ribozymes seems simpler, though.

    - Add more DNA the in vitro transcription reaction.  The promoter competes 
      with duplex RNA for T7 RNAP binding, so adding more promoter should 
      reduce the amount of 3' end extension.  That said, I already add about as 
      much DNA as I can to these reactions.

    - Optimize the sequence to disfavor 3' end extension, basically by removing 
      opportunities for duplex RNA formation near the 3' end.  The Y-tag 
      sequence is already pretty good in this regard, though, being mostly G.  

  - The 4-bp AAAA overhang on the mRNA is getting truncate somehow.  T7 RNAP 
    does occasionally stop a few nucleotides early [Gholamalipour2018]_, but 
    not nearly so often as it adds extra nucleotides, and not 50% of the time.  
    Possible solutions:

    - Using a cis-ribozyme would also solve this problem, since the end of the 
      transcript would occur much after the ribozyme.

  Either of the above hypotheses could be confirmed by RNA-seq, although that 
  sounds like more effort that this is worth.

Discussion
==========
- Adding excess linker does not increase yield.

- About half of the mRNA appears to be incapable of reacting with the linker.

- I could use 0.5x linker to reduce the amount of unligated linker without 
  decreasing the amount of ligated linker.
