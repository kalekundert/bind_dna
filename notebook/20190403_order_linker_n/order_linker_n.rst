**************
Order linker-N
**************

Kettner pointed out that there are companies that will synthesize complex 
oligos such as linker-N from [Naimudden2016]_, so there's no need to worry 
about making it myself.  He recommended `Midland CRC <oligos.com>`_, but I 
might be able to find others as well.

2019/04/03
==========
See attached quote from Midland CRC:

:download:`quotes/20190416_midland_63451.pdf`

.. update:: 2019/12/20

   I didn't order the complete linker-N.  I ordered just the branch with the 
   puromycin; they didn't synthesize the branch with the RT primer.  

   Well, now that I have to order more oligo anyways, I might as well find out 
   if things work better with Cy5 and 5' phosphorylation.

2020/02/20
==========
Since I have to reorder linker-N anyways, now is a good opportunity to consider 
if I want to make any changes to the oligo described by [Naimudden2016]_.  
Based on my experience up to this point, here are the things I'm thinking 
about:

Considerations
--------------
- Cy5 instead of FITC

  Most of the experiments I've been doing with linker-N involve gel 
  electrophoresis following by fluorescent imaging to see whether or not two 
  species (e.g. linker-N and mRNA, linker-N/mRNA and protein) are colocalized.  
  This requires imaging with two channels, which requires two fluorophores that 
  are well-separated, bright, and compatible with the laser scanners we have in 
  the lab.  So far, this has meant fluorophores excited by the 488 nm and 658 
  nm lasers.
  
  488 nm fluorophores:
  
  - FITC
  - GelRed/GelGreen
  - SYBR Green/Gold
    SYBR Green II (RNA)
  - SYPRO Orange
  - GFP

  658 nm fluorophores:

  - Cy5
  - SYPRO Ruby (some crosstalk, I think?)
  - Coomassie [Butt2013]_ (I haven't tried this yet).

  Having FITC in linker-N "uses up" the 488 nm channel, leaving me with fewer 
  options for staining proteins (although I just learned about Coomassie, that 
  might be promising) and no options for staining RNA.  In contrast, if 
  linker-N were labeled with Cy5, I would have good options for visualizing any 
  partner.  For example, :expt:`20200129_ligate_with_5_phosphorylation` uses a 
  Cy-5 labeled "pseudo-linker", which allows me to use the 488 nm channel to 
  visualize mRNA using GelRed.

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

- 5' phosphate

  The first step in cDNA display is to ligate linker-N to mRNA.  This requires 
  the 5' end of linker-N to be phosphorylated.  [Naimudden2016]_ perform the 
  ligation enzymatically, but there are a number of benefits to simply adding 
  the phosphate during synthesis:

  - This guarantees that every oligo is phosphorylated.  Enzymatic 
    phosphorylation is also efficient, but probably not perfect.  It's also 
    possible for the enzyme to go bad, making it an potential point of failure.

  - T4 PNK is sensitive to high-salt concentrations.  This causes headaches, 
    because the annealing step requires a high concentration of salt, and as a 
    result the annealing reaction either need to be dialyzed or significantly 
    diluted prior to phosphorylation, see 
    :expt:`20191219_anneal_linker_n_with_salt`.  Removing the phosphorylation 
    step eliminates this headache.

  - Synthetically-phosphorylated pseudo-linker seems to be ligated slightly 
    better than enzymatically-phosphorylated pseudo-linker, see 
    :expt:`20200129_ligate_with_5_phosphorylation`.  The difference is subtle, 
    and both are ligated well, but there is noticeably less of the 
    synthetically-phosphorylated linker left unligated.

  For normal oligos cost is a factor in deciding to do enzymatic vs. synthetic 
  phosphorylation, but for this oligo a single extra modification will be a 
  drop in the bucket.  It's also usually considered "worth it" to order 
  phosphorylated oligos for applications that require high efficiency, e.g.  
  library cloning.  I will be ultimately using linker-N to make libraries, and 
  I will care about efficiency, so I think phosphorylation is a no-brainer.

- Longer reverse-transcription (RT) primer arm

  I noticed that my pseudo-linker (o93) is ligated much more efficiently than 
  linker-N (o100), see 
  :expt:`20200213_ligate_linker_n_with_optimized_conditions`.  There are 
  several possible reasons for this (discussed in the linked experiment), but 
  one is that a longer RT primer arm is needed to keep the linker annealed.  
  The RT arm in linker-N is supposed to be `CCTTG`, but as noted above, this 
  arm is missing from the linker-N I ordered in April.  So I can't really say 
  if a longer arm than that is needed, but I can't really see how a longer arm 
  would hurt.

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
  
- Mass spectrometry (MS) quality control (QC)

  Another possible reason why my pseudo-linker might be ligated more 
  efficiently than linker-N is that linker-N wasn't synthesized correctly.  
  Note that I don't think that this is a particularly likely scenario, but it'd 
  be nice to have some data attesting that the oligo is what it should be.  

  [Naimudden2016]_ used MALDI-TOF to verify the identity of linker-N.  IDT, I 
  just learned, uses electrospray ionization-ion trap MS (ESI-IT) to verify 
  every oligo they synthesize.  They also have `an article 
  <https://www.idtdna.com/pages/education/decoded/article/esi-mass-spectrometry-why-we-use-it-for-oligonucleotide-quality-control>`__ 
  explaining that they find ESI-IT to be more accurate than MALDI-TOF for large 
  oligos.

  So I'm going to make sure I get MS QC from whichever company I end up 
  ordering from.  I think that will give me some peace of mind.  I'm also 
  asking all the companies to help me pick with MS method is most appropriate 
  for this oligo, since I don't really know anything about MS.

Inquiries
---------
Mostly I'm just keeping these in case I want to reuse the language to request 
more quotes.

2020/02/20 --- Inquiry to BEX:

   Below is the sequence of a branched oligo I would like to have synthesized.  
   Note that the "side chain" sequence is attached to the "main chain" via the 
   5-Me-dC brancher phosphoramidite:

   Main chain: 5′-(5' Phosphorylation)-CCCCCCCGCCGCCCCCCG-(5-Me-dC)-AAAAAAAAAAAAAAAAAA-(Spacer18)-(Spacer18)-(Spacer18)-(Cy5)-(Spacer18)-CC-(Puromycin)-3′
   Side chain: 5′-CCTTG-3′

   Please let me know if you can synthesize this oligo.  If you can, I'd like to get a quote for the smallest possible quantity I can order.  If the internal Cy5 is a problem, I'd like to get a quote with that phosphamidite replaced by (Fluorescein-dT) instead.

   Thank you,
   -Kale Kundert

2020/02/20 --- Inquiry to Midland CRC:

   Hi Sandy,

   Can you make me a quote for the smallest possible quantity of the following
   branched oligo?  Note that the "side chain" sequence is attached to the "main
   chain" via the 5-Me-dC brancher phosphoramidite, as I described previously:

   Main chain:
   5′-(Phosphate)-CCCCCCCGCCGCCCCCCG-(5-Me-dC)-AAAAAAAAAAAAAAAAAA-(Spacer18)-(Spacer18)-(Spacer18)-(Cy5)-(Spacer18)-CC-(Puromycin)-3′
   Side chain: 5′-CCTTG-3′

   Below are the Glen Research product numbers for the non-standard
   phosphoramidites referenced above:

   5-Me-dC: 10-1018
   Spacer18: 10-1918
   Cy5: 10-5915
   Puromycin: 20-4040

   I'm not sure what phosphoramidite is typically used to add a 5' phosphate, maybe
   10-1900, but hopefully you can let me know the best way to do that.  I'm also
   not sure if Cy5 can be incorporated internally as I have here.  If this can't be
   done, can you make me a quote with the Cy5 replaced by Fluorescein-dT (10-1056)
   instead?

   I also have some questions about your QC options.  Can you provide mass spec
   QC?  If so, what kinds of mass spec do you offer (e.g. MALDI, ESI, etc.) and
   which would be the most appropriate for an oligo of this size?  What degree of
   accuracy can I expect (e.g. mass within 10 Da, 1 Da, 0.1 Da)?  Please include
   the QC in the quote, as well.

   Thank you,
   -Kale Kundert
   
2020/02/21 --- Inquiry to `Biosynthesis <mailto:info@biosyn.com>`__:

   To whom it may concern,

   My name is Kale Kundert, and I'm a postdoc in the Church lab at Harvard.  I'd like to get a quote for the synthesis of the custom branched oligo given below.  Note that the "secondary chain" should be attached to the "primary chain" via the oxyhexyl group of the 5-Me-dC brancher phosphoramidite:

   Primary chain: 5'-(Phosphate)-CCCCCCCGCCGCCCCCCG-(5-Me-dC)-AAAAAAAAAAAAAAAAAA-(Spacer18)-(Spacer18)-(Spacer18)-(Cy5)-(Spacer18)-CC-(Puromycin)-3'
   Secondary chain: 5'-CCTTG-3'

   Below are the Glen Research catalog numbers for the non-standard phosphoramidites referenced above, in case the names I used aren't clear:

   5-Me-dC: 10-1018
   Spacer18: 10-1918
   Cy5: 10-5915
   Puromycin: 20-4040

   I'd also appreciate it if you could advise me on which of the mass spec QC method you offer (e.g. MALDI-TOF, LC-MS, ESI) would be most appropriate for an oligo of this size and complexity.

   Thank you,
   -Kale


Quotes
------
Midland: :download:`quotes/20200225_midland_63815.pdf`

.. datatable:: quotes/quotes.xlsx

- The given yield for Midland CRC is woth HPLC+PAGE purification, which is what 
  their scientists recommend.  Without PAGE, yield would be 5-10x greater, and 
  the price would be $125 less.


