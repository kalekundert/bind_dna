**************
o129, ribozyme
**************

:expt:`1` and :expt:`116` show that only about half of the mRNA appears to be 
capable of being ligated to the linker.  This maybe due to heterogeneity at the 
3' end of the mRNA.  There are several ways to eliminate that heterogeneity 
[Avis2012]_, but the simplest and most popular is to include a cis-acting 
ribozyme after the sequence of interest [Schurer2002]_ [Walker2003]_ 
[Akoopie2018]_.  The HDV ribozyme is efficient and puts no sequence 
requirements on the transcript.

Considerations
==============

Which ribozyme?
---------------
There are a couple different variants of the HDV ribozyme that have been 
described in the literature:

.. datatable:: hdv_ribozymes.xlsx

- [Schurer2002]_:

  - 67 nt
  - Probably the most widely used.
  - Least efficient if the nucleotide 5' of the cleavage site is an A, which it 
    is in my case.  But the effect isn't large (71% vs 79% cleavage).

- [Akoopie2018]_:

  - 56 nt
  - Reported to be better for some sequences and worse for others, but slightly 
    worse that [Schurer2002]_ on average.  Might be worth trying both.

I think I'll start with [Schurer2002]_.

Terminator
----------
Currently, I digest the plasmid encoding the mRNA to prevent RNAP from reading 
beyond the end of the Y-tag.  With a 3' ribozyme, this step may not be 
necessary.  Instead, I can add a T7 terminator after the ribozyme.  Termination 
won't be 100% effective, but it doesn't matter because every transcript will 
get processed to the right size anyways.  

While this would remove the need for the restriction digest, it might also 
introduce the need for a phenol-chloroform extraction (to remove RNase from 
miniprepped plasmid DNA).

.. update:: 2020/10/28

  In the interest of being able to clearly see if the ribozyme worked as 
  intended, I digested the plasmid with HaeII.  This should lead to the 
  creation of two fragments with defined length: the gene (829 bp) and the tail 
  (347, 336 bp for HDV67, HDV56).

Note that the whole terminator is transcribed, so I don't need to worry about 
including a spacer between the ribozyme and the terminator.  Briefly, 
rho-independent termination works like this [Peters2011]_:

- The terminator consists of a GC-rich hairpin followed by a U-rich tract.
- RNAP pauses on the U-rich tract.  I'm not totally sure why.
- This allows time for the GC-rich hairpin to form.
- The hairpin dislodges the RNAP.

Cloning
-------
Cloning in the Y-tag region is tricky, because there's a lot of GC- or AT-rich 
sequences.  I think the best bet is to order the HDV/T7-terminator cassette as 
a gBlock, then install it using Gibson assembly.  I decided to use XmnI to 
linearize the backbone.  This prevents the inserts from being specific to the 
affinity tag (e.g. His6, Strep, etc.) and avoids any difficult PCR steps.  It 
does keep the AAAA overhang created by XmnI, rather than the AAA overhang used 
by [Naimudden2016]_.  I don't think this overhang is a problem, though.  And in 
any case, it avoids changing two variables at once.

3' hydroxyl terminus
--------------------
As I had to discover the hard way, the 5' product of HDV ribozyme cleavage has 
a 2',3'-cyclic phosphate terminus.  Such a terminus cannot participate in 
ligation reactions, which require a 3' hydroxyl.

.. figure:: hdv_mechanism.png

  Mechanism of the HDV ribozyme.  Note that the 5' product is left with a 
  2',3'-cyclic phosphate terminus.  By Noxwei - Own work, CC BY-SA 3.0, 
  https://commons.wikimedia.org/w/index.php?curid=30227140

[Avis2012]_ actually addresses this issue very specifically:

   It is important to consider the chemical products of ribozyme cleavage as 
   both the HH Rz and HDV Rz yield a 5'-hydroxyl and a 2',3'-cyclic phosphate 
   (see Fig 2c).  The 5'-hydroxyl on a ribozyme-derived 5'-end is generally 
   more useful than its 5'-triphosphate counterpart.  Direct phosphorylation 
   of a ribozyme-derived 5'-hydroxyl is usually highly efficient since the 
   reaction is no longer reliant upon the successful removal of the 
   5'-triphosphate by alkaline phosphatase.  However, the 2',3'-cyclic 
   phosphate on a ribozyme-derived 3'-end is generally less useful than its 
   3'-hydroxyl counterpart.  It is important to be aware that RNA ligation or 
   3'-end radiolabeling (via ligation of 
   (5'-³²P)-cytidine-3',5'-bisphosphate) is not possible without first 
   removing the 2',3'-cyclic phosphate using the phosphatase activity of the 
   T4 PNK enzyme [Cameron1977]_, [Povirk1990]_.  It should also be noted that 
   an RNA with a 2',3'-cyclic phosphate carries additional charge and will 
   migrate approx 1 nt faster during electrophoresis than its 3'-hydroxyl 
   counterpart.  

[Schurer2002]_ references three specific protocols for recovering the 3' 
hydroxyl:

- Up to 50 pmol RNA were incubated with 6 U T4 polynucleotide kinase (New 
  England Biolabs) in 100 mM Tris–HCl pH 6.5, 100 mM magnesium acetate, 5 mM 
  β-mercaptoethanol in a final volume of 50 µl for 6 h at 37°C (14).

- Transcripts were incubated with 0.1 mM ATP, 100 mM imidazole–HCl pH 6.0, 10 
  mM MgCl2, 10 mM β-mercaptoethanol, 20 µg/ml BSA and 1 U T4 polynucleotide 
  kinase (New England Biolabs) per 100 pmol RNA in a final volume of 50 µl 
  for 6 h at 37°C (15).

- Up to 300 pmol tRNA were incubated in 100 mM morpholinoethanesulfonate-NaOH 
  pH 5.5, 10 mM MgCl2, 10 mM β-mercaptoethanol, 300 mM NaCl and 10 U T4 poly- 
  nucleotide kinase (New England Biolabs) in a final volume of 20 µl for 6 h 
  at 37°C [Povirk1990]_.

.. note::
  
   NEB T4 PNK buffer is:

   - 70 mM Tris-HCl
   - 10 mM MgCl2
   - 5 mM DTT
   - (pH 7.6 @ 25°C)

All three protocols are pretty similar: incubate the RNA with T4 PNK for 6h at 
37°C in a pretty standard buffer (no ATP).  I'm curious if I could just add the 
PNK to the transcription reaction and extend the reaction for 6h.


Preparation
===========
2021/06/21:

My Gibson assembly for p177/p178 didn't work on the first try.  I want to:

- Check that the XmnI digestion worked
- Try new Gibson mix.

.. figure:: 20210621_check_digestion_f120.svg

Observations:

- The digestion didn't go to completion.  This might be why I have a lot of 
  backbone colonies.

  How to deal with this:

  - Repeat digestion for longer time.
    
  - Gel purify: this would be hard; the bands aren't that separate.

  - Order new enzyme: I do use this enzyme for other things

  Let me start by just doing the digestion for longer, maybe 1h.

  NEB `troubleshooting guide <https://www.neb.com/tools-and-resources/troubleshooting-guides/restriction-enzyme-troubleshooting-guide>`__:

  - Use 5-10 U/µg DNA
  - Digest for 1-2h
  - Ensure that DNA < 25% total reaction volume

.. update:: 2021/06/23

  I got good results after making two changes:

  - I used my less-expired tube of 2x Gibson master mix.
  - I digested the backbone for 2h, and added extra enzyme after 1h30.

  I don't know which change was responsible.  It's possible that both my XmnI 
  and my original Gibson master mix are going bad.

Results
=======

−PNK --- 2020/10/28
-------------------
.. figure:: 20201028_ligate_with_ribozyme.svg

- The ligation reaction did not proceed because the 3' ends produces by HDV 
  ribozyme are not compatible with the ligation reaction.  See the "3'-hydroxyl 
  terminus" section above for a complete discussion of this.

- It's hard to say if the ribozyme reaction went to completion, although it 
  definitely proceeded to an appreciable extent.  The band representing the 
  cleaved ribozyme is clearly visible.  However, the band representing the 
  cleaved mRNA is so diffuse that it's hard to say for sure whether or not 
  there is a band representing the uncleaved transcript.

  I measured the intensities of the 829 bp bands relative to the ≈330 bp bands.  
  Accounting for the different lengths of the every RNA species, and assuming 
  that any excess intensity in the 829 bp band can be fully attributed to the 
  uncleaved transcript, it seems that only 60-70% of the transcripts were 
  cleaved.  This would be more convincing if the uncleaved bands was visible, 
  though.

  .. datatable:: 20201028_ligate_with_ribozyme.xlsx

- The HDV67 lanes have a third band at ≈200 bp.  I'm not sure what this band 
  is.

+PNK --- 2021/07/01
-------------------
.. figure:: 20210702_fix_ends_f123_f124.svg

Observations:

- The bands didn't run as far as they usually do.  Compare the ladder from this 
  gel to the ladder in the experiment above: the 150bp band should be about 50% 
  down the gel,  but is only 25% here.

  I think the problem was that the current setting on the gel was much too low, 
  so the voltage was lower than usually because the current was capped.  Ugh...

- I think the HDV56 ribozyme might cleave itself more efficiently, but it's 
  hard to be sure.  I only see the cleaved ribozyme band for the HDV56 
  ribozyme, and although I think I see the cleaved mRNA bands for both 
  ribozymes, it appears stronger relative to the uncleaved band for HDV56.

- I should run a gel of the purified mRNA.  Not sure if I just do that as part 
  of the ligation reaction.  Even though f124 seems to be cleaved more 
  efficiently, I recovered more f123.  That may indicate that my gel 
  purification wasn't very good.

- I need to use razors to cut out these bands.  Those plastic tube things get 
  too much.


.. protocol:: 20210707_ligate_with_ribozyme.pdf 20210707_ligate_with_ribozyme.txt
.. figure:: 20210707_ligate_with_ribozyme.svg

Observations:

- The gel purification described above was not very effective.  A lot of 
  uncleaved mRNA and cleaved ribozyme was present in both purified samples.

- The ribozyme reaction apparently proceeded to completion during the 
  annealing/ligation reactions.  This is surprising to me because it didn't go 
  to completion after 4h transcription and 6h KLD reactions.  Some steps that 
  could be responsible:

  - The 90°C incubation and subsequent ramp-down.
  - The T4 RNA ligase buffer (50 mM Tris, 10 mM MgCl₂, 1 mM DTT, pH 7.5)

  I think the 90°C incubation is the more likely culprit.  The T4 RNA ligase 
  buffer is almost identical to the T4 PNK buffer used in the dephosphorylation 
  reaction (70 mM Tris, 10 mM MgCl₂, 5 mM DTT, pH 7.6).

  This is interesting because one of the concerns with the ribozyme approach is 
  that the cleavage step itself is relatively low efficiency.  This seems to be 
  a solution to that problem.

- The ligation efficiency seems low.  It's hard to say exactly what the 
  ligation efficiency is, because there is clearly a lot of freshly-cleaved 
  (and therefore non-PNK-treated mRNA) in the reaction.  But just comparing the 
  amount of ligated mRNA to the amount of KLD-treated mRNA in the −linker 
  control, it's clear that only a small fraction of the mRNA was ligated.

  Poor efficiency could be attributed to T4 PNK or T4 RNA ligase.  Both should 
  work well, though.  [Schurer2002]_ claims that dephosphorylation goes to 
  completion, and I've seen in the past that ligation is efficient (for the 
  subset of mRNA that is reactive).

- More of the HDV56 mRNA is ligated, compared to the HDV67 mRNA.  I think this 
  is just because the HDV56 purification recovered more of the cleaved and 
  dephosphorylated starting material.

- There are some bands I can't identify:
  
  - The two bands above the 411 nt cleaved mRNA.  Both disappear in the +linker 
    lanes, suggesting that both are uncleaved mRNA+HDV.  The upper one is close 
    to the expected MW for this band.  The lower one is two small

  - The 100 nt bands in the +linker conditions.  This is probably cleaved 
    ribozyme, given that it's not present in the −linker condition.

  I think it's likely that my DNA template didn't amplify very cleanly.  I 
  never ran a gel of the template, so I don't know what it looks like.  But I 
  seem to have two dominant template lengths.  Both would contain full-length 
  ribozyme, because I included some buffer to make the cleaved HDV band a 
  certain size.

Conclusions:

- If I can get cleavage to go to completion, it would make sense to add a 
  restriction site just after the ribozyme.

  - Making the ribozyme as short as possible would make gel purification 
    easier.
  - Linearization by restriction digest is less error-prone than linearization 
    by PCR.

- I'm inclined to try this again.  I think I can do a much better job making 
  the starting material—especially if the temperature trick works—and its 
  possible that the efficiency will look much better with clean starting 
  material.
