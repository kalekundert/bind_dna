**************
Order linker-N
**************

Kettner pointed out that there are companies that will synthesize complex 
oligos such as linker-N from [Naimudden2016]_, so there's no need to worry 
about making it myself.  He recommended `Midland CRC <oligos.com>`_, but I 
might be able to find others as well.

Considerations
==============

Cy5 instead of FITC --- 2020/02/20
----------------------------------
Most of the experiments I've been doing with linker-N involve gel 
electrophoresis following by fluorescent imaging to see whether or not two 
species (e.g. linker-N and mRNA, linker-N/mRNA and protein) are colocalized.  
This requires imaging with two channels, which requires two fluorophores that 
are well-separated, bright, and compatible with the laser scanners we have in 
the lab.  So far, this has meant fluorophores excited by the 488 nm and 658 nm 
lasers.

488 nm fluorophores:

- FITC
- GelRed/GelGreen
- SYBR Green/Gold
- SYBR Green II (RNA)
- SYPRO Orange
- GFP

658 nm fluorophores:

- Cy5
- SYPRO Ruby (some crosstalk, I think?)
- Coomassie [Butt2013]_ (seems to quench 488 nm fluorophores)

Having FITC in linker-N "uses up" the 488 nm channel, leaving me with fewer 
options for staining proteins (although I just learned about Coomassie, that 
might be promising) and no options for staining RNA.  In contrast, if 
linker-N were labeled with Cy5, I would have good options for visualizing any 
partner.  For example, :expt:`16` uses a Cy-5 labeled "pseudo-linker", which 
allows me to use the 488 nm channel to visualize mRNA using GelRed.

I'm not totally sure that Midland CRC will be able to install Cy5 internally, 
though.  I know they can use any phosphoramidite sold by Glen Research, and 
the Cy5 phosphoramidite is catalog number 10-5915.  This phosphoramidite does 
have a "5'" alcohol, but it's protected by MMT instead of DBT.  Glen Research
notes that "cyanine dyes are normally added once at the 5'-terminus", but 
clearly it is possible to deprotect the Cy5 phosphoramidite and expose a "5'" 
nucleophile, even if the deprotection conditions might be different than 
usual (although MMT and DBT are very similar chemically).  It's also worth 
noting that IDT sells oligos with internal Cy5, although most other vendors 
I've seen only explicitly offer 5' Cy5.

In contrast, FITC is incorporated into linker-N via a dT with FITC attached 
to the base (10-1056).  This can be used internally with no issues at all.  I 
believe it doesn't even interfere with base pairing, although that's not 
relevant to linker-N.

Assuming that internal Cy5 is an option, I'll get it.  Otherwise I'll stick 
with FITC and just make it work with the two-channel images.  It should be 
less of a problem once I finish debugging the ligation reaction and move on 
to the expression/puromycin reaction.

.. update:: 2020/03/09

   It seems that all vendors can install internal Cy5 modifications, so 
   that's what I'm going to order.

5' phosphate --- 2020/02/20
---------------------------
The first step in cDNA display is to ligate linker-N to mRNA.  This requires 
the 5' end of linker-N to be phosphorylated.  [Naimudden2016]_ perform the 
ligation enzymatically, but there are a number of benefits to simply adding the 
phosphate during synthesis:

- This guarantees that every oligo is phosphorylated.  Enzymatic 
  phosphorylation is also efficient, but probably not perfect.  It's also 
  possible for the enzyme to go bad, making it an potential point of failure.

- T4 PNK is sensitive to high-salt concentrations.  This causes headaches, 
  because the annealing step requires a high concentration of salt, and as a 
  result the annealing reaction either need to be dialyzed or significantly 
  diluted prior to phosphorylation, see :expt:`12`.  Removing the 
  phosphorylation step eliminates this headache.

- Synthetically-phosphorylated pseudo-linker seems to be ligated slightly 
  better than enzymatically-phosphorylated pseudo-linker, see :expt:`16`.  The 
  difference is subtle, and both are ligated well, but there is noticeably less 
  of the synthetically-phosphorylated linker left unligated.

One possible reason not to add the 5' phosphorylation during synthesis is that 
it may interfere with synthesis of the RT arm.  [Naimudden2016]_ protects the 
5'-OH of the ligation arm with an acetyl group during synthesis of the RT arm.  
This makes sense; the 5'-OH would definitely interfere with synthesis.  But I 
don't know if the phosphate group needs protection, or what groups could 
protect it.  It's possible no one's tried to make anything like this before.  

For normal oligos cost is a factor in deciding to do enzymatic vs. synthetic 
phosphorylation, but for this oligo a single extra modification will be a 
drop in the bucket.  It's also usually considered "worth it" to order 
phosphorylated oligos for applications that require high efficiency, e.g.  
library cloning.  I will be ultimately using linker-N to make libraries, and I 
will care about efficiency, so I think phosphorylation is a no-brainer 
(assuming it doesn't interfere with RT-arm synthesis).

RT primer 3'-OH --- 2020/03/09
------------------------------
The RT primer arm needs to end with a 3'-OH in order to be able to prime RT.  
However, normal oligo synthesis goes in the 3'→5' direction, which would result 
in the arm ending with a 5'-OH.  (To put it another way, normal oligo synthesis 
can only have one 3'-OH.  Each branch only creates another 5'-OH.)  To account 
for this, the RT arm needs to be synthesized using 5' phosphoramidites rather 
than the usual 3' phosphoramidites.  These are available from `Glen research 
<https://www.glenresearch.com/applications/specialized-dna-and-rna-synthesis/5-3-synthesis-phosphoramidites-and-supports.html>`__:

- dA-5'-CE Phosphoramidite (10-0001)
- dC-5'-CE Phosphoramidite (10-0101)
- dT-5'-CE Phosphoramidite (10-0301)
- dmf-dG-5'-CE Phosphoramidite (10-9201)

  .. note::
   
     dG-5' has a different protecting group (dmf) than "normal" dG (iBu).  I 
     was initially concerned because I thought that dmf was some kind of 
     modification.  This `appears 
     <https://www.glenresearch.com/reports/gr9-12>`__ not to be the case, and I 
     don't think there's any reason not to use dmf-dG.

For comparison, here are the product numbers for the normal 3'→5' oligos:

- dA-CE Phosphoramidite (Glen 10-1000)
- dC-CE Phosphoramidite (Glen 10-1010)
- dG-CE Phosphoramidite (Glen 10-1020)
- dT-CE Phosphoramidite (Glen 10-1030)

Note that synthesis needs to begin with 3'-puromycin, since the puromycin 
monomer is only available attached to a solid support.  That's why the RT and 
ligation arms can't be contiguous.

Longer RT primer arm --- 2020/02/20
-----------------------------------
I noticed that my pseudo-linker (o93) is ligated much more efficiently than 
linker-N (o100), see :expt:`17`.  There are several possible reasons for this 
(discussed in the linked experiment), but one is that a longer RT primer arm is 
needed to keep the linker annealed.  The RT arm in linker-N is supposed to be 
`CCTTG`, but as noted above, this arm is missing from the linker-N I ordered in 
April.  So I can't really say if a longer arm than that is needed, but I can't 
really see how a longer arm would hurt.

It's worth nothing that the 5 nt RT arm is only 1 nt shorter than the 
random-hexamers that are often used to prime RT reactions.  Of course, the RT 
arm is also held in place by 17 nt of complementarity on the other side of 
the puromycin arm, so it should be well-anchored.

The question is whether the puromycin arm, along with the 5-Me-dC brancher, 
will interfere with RT binding/function.  I tried to see if I could find a 
structure of the MMLV RT (which is what [Naimudden2016]_ use) in complex with 
DNA/RNA, so get a sense for how long the puromycin arm would need to be to 
*not* be in the way.  There are some structures [Cote2008]_, but it seems 
necessary to piece together information from a number of different structures 
and experiments to say anything about how MMLV RT binds DNA, and I don't 
think I have to domain knowledge necessary to do that.

I think I'm going to leave the arm as-is for now.  I'm pretty resistant to 
making any changes to the sequence of linker-N, because I know the sequence 
published by [Naimudden2016]_ should work.  I just can't really be sure that 
any change will be an improvement, even if it makes sense.  And in this case 
the justification for a change isn't very strong.

Mass spectrometry --- 2020/02/20
--------------------------------
Another possible reason why my pseudo-linker might be ligated more efficiently 
than linker-N is that linker-N wasn't synthesized correctly.  Note that I don't 
think that this is a particularly likely scenario, but it'd be nice to have 
some data attesting that the oligo is what it should be.  

[Naimudden2016]_ used MALDI-TOF to verify the identity of linker-N.  IDT, I 
just learned, uses electrospray ionization-ion trap MS (ESI-IT) to verify 
every oligo they synthesize.  They also have `an article 
<https://www.idtdna.com/pages/education/decoded/article/esi-mass-spectrometry-why-we-use-it-for-oligonucleotide-quality-control>`__ 
explaining that they find ESI-IT to be more accurate than MALDI-TOF for large 
oligos.

So I'm going to make sure I get MS QC from whichever company I end up 
ordering from.  I think that will give me some peace of mind.  I'm also 
asking all the companies to help me pick with MS method is most appropriate for 
this oligo, since I don't really know anything about MS.  I'll be fine with 
either TOF or ESI, though.

Purification --- 2020/03/09
---------------------------
There are 3 commonly-used methods for purifying complex oligos.  The pros and 
cons of these methods are described in `this article 
<https://www.sigmaaldrich.com/technical-documents/articles/biology/best-purification.html>`__, 
and summarized here:

- Reverse-phase HPLC (RP-HPLC): RP-HPLC separates oligos based on differences 
  in hydrophobicity.  This is an effective way to purify full-length oligos, 
  because only full length oligos will have a 5'-DMT protecting group (very 
  hydrophobic) at the end of synthesis.  (The 5'-DMT is removed after, or 
  during, purification.)  RP-HPLC is also effective for purifying oligos with 
  hydrophobic modifications (e.g.  dyes), because both the modification and the 
  5' DMT will affect mobility.  That said, linker-N has so many hydrophobic 
  spacers that I might end up selecting only for that; the 5'-DMT and iCy5 
  might become relatively insignificant.
  
- Ion exchange HPLC (IE-HPLC): IE-HPLC separates oligos based on charge (i.e.  
  in the phosphate backbone).  This method is most effective for oligos shorter 
  than 40 nt.  For longer oligos, truncated products will have enough charge to 
  be not-easily-distinguishable from full-length product.  IE-HPLC is also 
  especially effective for oligos with significant secondary structure, because 
  it is compatible with mobile phases that disrupt base-pairing.  I don't think 
  this method is appropriate for linker-N.

- PAGE: PAGE separates oligos based on size and charge.  This is generally 
  considered the best purification method, but it also results in the lowest 
  yields.  Because size directly affects migration through the gel, PAGE 
  purification works well regardless of how hydrophobic/charged the oligo is.  

  When I've run o100 on PAGE, I've seen multiple bands.  The same bands are 
  evident in [Naimudden2016]_.  More specifically, there are a number of bands 
  around the expected size, then one band about double that size.  I've assumed 
  that these were all different ways for the linkers to anneal or otherwise 
  stick together, but now that I think about it, that explanation can really 
  only apply to the highest band.  The rest are probably synthesis errors.  
  Even the highest band could be an error, e.g. shortmers consisting of little 
  more that the 3' spacers may not migrate very fast because they wouldn't be 
  very charged.  So I have reason to suspect that PAGE purification would give 
  me substantially more pure linker.  I'm also curious what the PAGE gel of 
  PAGE-purified linker-N would look like.

Probably HPLC purification would give better yield, while PAGE purification 
would give better purity.  To choose between these methods, I need to decide 
how important purity is to me.  In general, I would probably say that purity is 
not that important.  Although errors could reduce the efficiency of the display 
reaction, they would not otherwise affect the results of my binding assay.  
Still, it's worth considering the specific errors that could occur and what 
effect of efficiency they might have:

- Most of the errors will be shortmers, meaning that synthesis was prematurely 
  terminated and the 5'-end of the oligo (3'-end in the reverse-synthesized RT 
  arm) is missing.  These errors will mostly preclude successful display.  
  Either the 5' phosphate will be missing, in which case the ligation will 
  fail, or part of the RT arm will be missing.  The RT arm is already short, so 
  I would expect any deletions to have a significant effect of priming.  Also, 
  missing the RT arm seems to affect ligation, so oligos with RT arm deletions 
  may not even be ligated well.

- Indels in the poly-A region wouldn't have any effect, because the length of 
  the spacer isn't critically important, and I'm not planning to do a poly-A 
  capture anyways.

- Indels in the spacers may not have much effect, because [Naimudden2016]_ 
  showed that the mRNA/protein fusion can form without any spacers at all 
  (albeit at lower efficiency).  
  
- Indels in the ligation arm may also be tolerated, because the arm is pretty 
  long.  That said, I already have evidence that o100 doesn't anneal/ligate 
  particularly well, so it may be that even a single mistake in this arm would 
  significantly impair ligation.

With the above in mind, here are some reasons for choosing one kind of 
purification over another:

- All told, I think it's safe to assume that most synthesis errors will be 
  non-functional.  This is an argument for doing a more rigorous purification, 
  since then my concentrations will be more reliable (e.g. I won't have to 
  worry that 50 pmol of linker only contains 10 pmol of *functional* linker) 
  and I won't be losing as much material (since more of the material I'm losing 
  wouldn't have worked anyways).

- I'm still trying to debug things, and the less things that can go wrong, the 
  easier that will be.  This is another argument for doing a more rigorous 
  purification, at least for now.  It might make sense to try a less rigorous 
  purification when I'm scaling up and need more material.

- [Naimudden2016]_ used RP-HPLC purification, so I know that should work.

- Maybe there's an argument for doing two purification (e.g. HPLC+PAGE), but 
  given that [Naimudden2016]_ only did one, I don't think I'll need to go this 
  far.

- It's a little hard to say whether PAGE of RP-HPLC would give the better 
  purification.  On one hand, IDT `recommends 
  <https://www.idtdna.com/pages/support/faqs/when-should-i-choose-page-vs-hplc-purification->`__ 
  RP-HPLC for modified oligos (because modifications often increase 
  hydrophobicity).  On the other hand, RP-HPLC may be less effective for an 
  oligo with this many hydrophobic modifications.  My PAGE gels suggest that 
  HPLC purification doesn't yield particularly pure oligo.

I think it would be best to get PAGE purification for now.  I should revisit 
this if I decide to scale up and need more linker.

Orders
======

2019/04/03
----------
See attached quote from Midland CRC:

:download:`quotes/20190416_midland_63451.pdf`

.. update:: 2019/12/20

   I didn't order the complete linker-N.  I ordered just the branch with the 
   puromycin; they didn't synthesize the branch with the RT primer.  

   Well, now that I have to order more oligo anyways, I might as well find out 
   if things work better with Cy5 and 5' phosphorylation.

2020/02/20
----------
Since I have to reorder linker-N anyways, now is a good opportunity to consider 
if I want to make any changes to the oligo described by [Naimudden2016]_.  See 
the Considerations_ section for more discussion.

Quotes:

.. datatable:: quotes/quotes.xlsx

- The yields are estimates, so should be taken with a grain of salt.  In truth, 
  I think I'd probably get about the same yield from any of these companies, 
  since it seems like they're all starting with about 250 nmol.  The 
  purification strategies are a little different (e.g. BEX proposes to do the 
  least purification, so they may give the best yield), but I'm not really sure 
  whether yield or purity is more important to me.

  .. update:: 2020/03/09

      I now think that PAGE purification is what I want.

- I might not even use all of the oligo, so it's more important to pick the 
  lowest price than the best value.  Especially considering that the value 
  depends on the yield, which may not be estimated accurately.


