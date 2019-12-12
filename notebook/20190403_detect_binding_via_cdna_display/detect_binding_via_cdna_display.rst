*******************************
Detect binding via cDNA display
*******************************

My first experiment will be to establish whether or not ligation can be used to 
detect binding between a protein and its DNA target.  This experiment will also 
serve as a test of the cDNA-display technology, which may be difficult to get 
working in its own right.  cDNA display is mRNA display (i.e. mRNA with a 
puromycin linker) with a reverse translation step.

- Advantages:

   - Most minimal linker between protein and gene.

   - Covalent linkage.

- Disadvantages:

   - Puromycin linker is either expensive to purchase or difficult to 
     synthesize.

   - Puromycin can induce premature chain termination:
      
      - This can happen either if:

         - The puromycin linked to the mRNA being transcribed gets in the 
           ribosome too early.

         - A puromycin from another mRNA gets in the ribosome.

      - Can purify for full length product, but this reduces yield.

      - Longer proteins are more likely to be prematurely terminated.  This 
        places a practical limit of 300 aa on the size of proteins that can be 
        expressed.

- References:

   - [Naimudden2016]_: Forked linker; latest in series of optimizations.

   - [Kurz2001]_: Linker attached to mRNA via psoralen crosslinking.


Considerations
==============

Reagents
--------
This is a consolidated and categorized list of all the reagents I will need to 
order for this experiment:

- General supplies

   - Eppendorf Safe-Lock Tubes, 1.5 mL, Eppendorf Quality™, colorless, 500 
     tubes (Eppendorf 022363204)

- Cloning

   - Esp3I (R0734S)

   - SapI (R0569S)

   - T4 DNA Ligase (NEB M0202)

   - MultiShot™ StripWell Mach1™ T1 Phage-Resistant Chemically Competent E.  
     coli (Invitrogen C869601)

      - These are really expensive, but very fast growing.  We'll see how I 
        feel about burning money, but my intention is just to make my own from 
        one of these aliquots.

   - QIAprep Spin Miniprep Kit (Qiagen 27104)

- To transcribe the mRNA:

   - HiScribe T7 Quick High Yield RNA Synthesis Kit (NEB E2040)

   - RNaseZap (AM 9780)

   - DNase/RNase-free water (Zymo W1001-10)

   - RNA Clean & Concentrator-5 (Zymo R1013, $130)

      - Decided against: RNeasy MinElute Cleanup Kit (Qiagen 74204, $360)

        More expensive.  But I wasn't thinking that it might be compatible with 
        the vacuum thing, that might be nice.

   - 6% TBE/Urea gels (Invitrogen EC68652BOX)

   - PAGE GelRed (Biotium 41008-T)

- To synthesize the cDNA linker:

   - https://secure.ids5.com/oligos/default.asp ($1500)

- To ligate the linker to the mRNA:

   - [Naimudden2016]_:

      - T4 polynucleotide kinase (NEB M0201)

      - T4 RNA ligase (Takara 2050B, $679)

   - [Kurz2001]_:

      - T4 DNA ligase

      - 450 W medium pressure immersion lamp (ACE Glass), equipped with a Pyrex 
        absorption sleeve in a Quartz immersion well.

- To express the protein:

   - PURExpress In Vitro Protein Synthesis Kit (NEB E6800)

   - PURExpress® Δ RF123 Kit (NEB E6850S)

   - Ni-NTA Magnetic Agarose Beads (Qiagen 36111)

   - Reverse transcriptase

      - Probably want an MMLV RNase H⁻ enzyme, for least RNase activity.

      - SMART® MMLV Reverse Transcriptase (Takara 639524) [Naimudden2016]_

   - RNase H (Ambion AM2293)

- To make double-stranded cDNA:

   - Polymerase

   - EcoRV-HF (NEB R3195)

   - BstNI (NEB R0168S)

   - NdeI (NEB R0111S)

   - BamHI-HF (NEB R3136S)

   - 5' protein barcode:

      - Phosphorylated primer

   - 3' protein barcode:

      - Unphosphorylated primer

- To perform the ligation assay:

   - PCR primers compatible with target and cDNA

- To quantify ligation by qPCR:

   - qPCR master mix

   - Plates

   - Seals

   - Reference amplicon (same length, different primers, known concentration)

   - Validated primers

      - Order expected products in advance.

      - Validate primers while waiting for other things.

Clone the protein
-----------------
There are a number of decisions to make regarding how to make the protein 
construct:

5' vs 3'  barcode
~~~~~~~~~~~~~~~~~
The barcode identifying the protein can either be placed before (5') or after 
(3') of the gene encoding the actual protein.  There are a lot of pros and 
cons, so it would probably be prudent to try both approaches.

3' barcode:

The advantage of the 3' barcode is that it would let me cleave off all of the 
cDNA but the barcode, which would ameliorate the problems discussed in the 5' 
barcode section below.  The problem is that it might be difficult to avoid the 
3' barcode from being translated.  Simply letting the barcode be translated 
would be bad:

- The barcode could affect the function of the protein.  For example, an 
  especially hydrophobic tag could destabilize the protein fold, or a 
  negatively charged tag could repel DNA.

- Each protein would (of course) have a different tag, so any effect the tags 
  have wouldn't be consistent between proteins.

- I could possibly test the effect of the barcodes on a control protein, but 
  even that wouldn't really be informative.  It's very possible that the effect 
  of the tags would depend on the specific protein it's linked to.

One way to avoid translation of the barcode would be to add a stop codon and 
translate in the absence of release factors:

- NEB has a PURExpress kit lacking release factors:

   - PURExpress® Δ RF123 Kit (NEB E6850)

- Stop codons are not recognized by tRNAs, but by "release factors" (which 
  are proteins).  So presumably, if the ribosome encountered a stop codon in 
  the absence of release factors, the A-site would just sit empty and 
  puromycin (if it were close enough) would be able to bind.
  
- The cDNA-display linkers have been optimized so that the puromycin is 
  correctly positioned to enter the A-site as the ribosome stalls where the 
  mRNA is ligated to the linker.  By causing the ribosome to stall earlier, 
  it might be necessary to repeat this optimization, which would probably be 
  difficult and time-consuming.

- If the ribosomes read through the stop codon at a significant enough rate, I 
  would have to do something about that.  Maybe add a pulldown or cleavage tag 
  after the stop codon, so I can remove proteins with barcodes expressed.
  
Another way to accomplish the goal of having the protein labeled only with a 
barcode would be to attach the DNA using an emulsion-based technique, e.g.  
[Yonezawa2003]_.  I haven't looked into this carefully, but basically since 
things are encapsulated in droplets, you have a lot more flexibility in how you 
digest things.
  
Another possible problem is that with just the barcode, the cDNA might be short 
enough that it would have trouble ligating with bound DNA:

- This could make ligation efficiency dependent on the orientation of the 
  protein binding domain relative to its C-terminus (where the cDNA would be 
  attached).  This would definitely not be desirable.

- This might be mitigated by the puromycin linker, which anyway needs to 
  contain a region long- and flexible-enough to reach the A-site from wherever 
  the mRNA is.

I had the mistaken idea that I could put a TEV site (or similar, e.g. IMPACT) 
before the barcode, then just cleave the barcode off after translation.  The 
problem here is that the mRNA is attached to the C-terminus, so if I cleave off 
a C-terminal tag, I'll lose the mRNA.  I could imagine putting the barcode in 
the middle of an intein.  In this way, the intein would cut itself and the 
barcode out, leaving the protein attached to its mRNA.  But the barcode again 
could have unpredictable effects on intein function, and this would be hard to 
control for.  Also, according to [Gu2014]_, ribosome display is limited to 
proteins of about 300 aa or less.  The IMPACT intein is 198 aa, which leaves 
only about 100 aa for my DNA binding domain (Zif268 is about 90 aa).

5' barcode:

The advantage of the 5' barcode is that the barcode is never translated, as 
discussed above.  The disadvantage is that the entire cDNA will be present for 
the binding reaction.  This could cause the following problems, which can be 
mitigated in various ways:

- The cDNA has a very high effective concentration relative to the protein it's 
  displaying.  As that protein is a DNA-binding protein, it might be difficult 
  for the targets to out-compete the cDNA itself.

   - However, this might also improve my signal-to-noise ratio by filtering out 
     weak binding events.

   - If this is a problem, it could be mitigated by adding more target DNA.  
     The target DNA will never match the local concentration of the cDNA, but 
     more target might help.

   - Removing the coding DNA (i.e. the advantage of the 3' barcode) might not 
     actually solve this problem, although it certainly wouldn't hurt.

- Proteins may find targets to bind in each others cDNAs.  This could result in 
  targets being ligated not to the cDNA of the proteins binding them, but to 
  the cDNAs of other proteins binding the cDNA of the protein binding the 
  target.

   - Keep the proteins dilute relative to the targets.

   - If using proteins with partially known targets, reverse translate the 
     proteins such that the cDNA doesn't contain any potential binding sites.

   - Do a control where a known target site is explicitly included in the mRNA, 
     and see how much cross-ligation occurs.  This could be part of a series of 
     experiments to determine a good protein:target:ligase ratio.

   - Just don't worry about it.  Most of the cDNA will be the same for most of the 
     proteins, so the effect of a protein that targets the cDNA will most likely 
     be too diffuse to matter (except for not getting a good signal for that 
     protein).

Restriction digest
~~~~~~~~~~~~~~~~~~
If the barcode is on the 5' end of the mRNA, there's no specific need to digest 
the cDNA.  However, wherever the barcode is, there are some advantages to 
digesting the cDNA:

- T7 polymerase can append variable numbers of G's to the beginning of 
  transcripts [Imburgio2000]_.  This variability might make it hard to 
  interpret 5' barcodes.  Adding a restriction site (or really any fixed non-G 
  sequence) would make interpreting the barcodes more reliable.
  
- A digest would naturally phosphorylate the end of the cDNA, which is 
  necessary for ligation.

- Could leave overhangs, which may be important depending on the ligation 
  strategy.  See the Ligation_ subsection for a more detailed discussion.

For my first construct, I decided to include a panel of restriction sites to 
allow me to experiment with different sticky end lengths:

- I decided to just use the same enzymes as [Bauer2017]_; it seems like a good 
  mix of robust enzymes with different overhangs:

   - EcoRV: Blunt-end, A/T
   - NruI: Blunt-end, G/C
   - BstNI: 1-bp overhang, A/T
   - NdeI: 2-bp overhang, AT/TA
   - BamHI: 4-bp overhang, GATC/CTAG

Barcode sequence
~~~~~~~~~~~~~~~~
Obviously I don't really need to barcode the protein in this assay, since 
there's only one protein.  But in order to test the assay most realistically, I 
want to use a barcode of the correct sequence and length.

Not surprisingly, there's plenty of literature on how to construct good 
barcodes.  A good barcode should:

- Be able to correct for a small number of sequencing/synthesis errors.
- Avoid long runs of homopolymers.
- Have relatively balanced GC content.
- Avoid sequences known to induce sequencing errors.

I decided to follow [Hawkins2018]_, because it seemed like a thoughtful, 
modern, and applicable approach.  The code to generate barcodes is available 
from `github <https://github.com/finkelsteinlab/freebarcodes.git>`_, but for 
now, I'm just using Barcodes17-2 (i.e. 17 nt long, capable of correcting 2 
errors) in the spreadsheet from the supplement.  This set includes 23025 
barcodes.  It's possible that I'll want more than that, but this should be a 
reasonable starting point and it doesn't require me to run any code.

The first barcode in this set is ``AACAACAACAACAACCG``, so that's what I'll use 
for this experiment.

Clone the target
----------------
There are several things to consider in the design of the target DNA molecule:

PCR primers
~~~~~~~~~~~
Ligation with the cDNA needs to create an amplicon that includes the target 
sequence.  There are two ways to approach this:

- One external primer:

   - ``[fwd]—[target]``

   - Since I couldn't amplify the target by PCR, I'd have to prepare in another 
     way:
     
      - Amplify it in a plasmid, then excise it using a restriction enzyme.

         - Actually, I might prefer to do things this way anyway.  PCR is 
           messy!  A digest would require a purification, but yield isn't 
           really a concern here (these reactions will be so small).

         - I want my targets to be dephosphorylated, so if I excise them via 
           restriction digest, I'll have to add a dephosphorylase.

         - I could design oligos such that each contains ~5 of these target 
           motifs, separated by restriction sites.  I couldn't do this if I were 
           going to amplify with PCR, because I'd get mixed products.

         - It's possible that restriction cloning would cause some plasmids to 
           harbor multiple targets, but they would all be freed by the digest.

      - Anneal oligos:

         - Not sure how efficient this is.  In other words, how much ssDNA 
           would I end up with?

         - Wouldn't be compatible with degenerate nucleotides.

   - I would need to use sticky ends to ensure that cDNA is ligated 3' of the 
     target.  This could increase the chance of spurious ligations.  I'd need 
     to test different overhang lengths to find an optimal combination of high 
     directional ligation and low inter-molecular ligation.

   - I'd be concerned that there isn't enough room for the ligase to operate 
     after the target, if the target is being bound by a library member.  I 
     might have to experiment with adding extra sequences 3' of the target.

   - Smallest sequence.

- Two external primers:
   
   - ``[fwd]—[target]—[rev]``

   - Because both primers point towards the target, the target will be 
     amplified no matter which primer is used (in conjunction with a primer 
     in the cDNA) to amplify successful ligations.

   - With blunt-ended ligation, I would only be able to amplify half of any 
     ligations.

      - Because the ligation can occur in either orientation with equal 
        probability, half of the products would need to be amplified with 
        the forward primer while the other half would need the reverse 
        primer. 

      - I can't add both primers to the same reaction, otherwise I'll just 
        amplify the targets.  It'd be a mess.

      - So I'd have to just pick a primer and accept that I'm losing half my 
        data.  This could limit my ability to detect weak interactions.

      - Using sticky ends (see above) would alleviate this problem, but might 
        also increase the chance of spurious ligations.  
        
- One internal primer:

   - ``[target]—[primer]—[target]``

   - Because there have to be at least two copies of the target, I couldn't 
     make targets using degenerate nucleotides.  I don't expect this to be 
     particularly limiting, though.
   
   - This layout provides ways to amplify blunt-end ligation products:

      - If the primer is palindromic, a single primer could be to used to 
        amplify "towards" the cDNA cassette regardless of the orientation of 
        the ligation.

         - PCR with palindromic primers may be inefficient, because the primers 
           can anneal with each other.  However, this wouldn't be as bad as 
           canonical "primer dimers", which are primers that only anneal 
           partially.  The overhangs from these primer dimer are filled in by 
           the polymerase and amplified, ultimately using up a lot of the PCR 
           reagents and creating a small product that needs to be purified 
           away.           

         - It would not be possible to use PCR to add sequencing adaptors.

         - A palindromic sequence with propensity to form hairpins may also 
           cause problems for sequencing.

      - If the primer isn't palindromic, I could add both the primer and its 
        reverse complement to the same reaction.
         
         - The PCR may still be inefficient as described above, because the 
           primers could still anneal with each other.

         - I again couldn't use PCR to add sequencing adaptors, because that 
           would result in the creation of primer dimers.

         - However, the final product wouldn't have any palindromic sequences
        
   - Would need to amplify via plasmid (see above).

   - Better target-to-primer ratio, compared to above.

- Two internal primers:

   - ``[target]—[rev]—[fwd]—[target]``

   - Could amplify all blunt-end ligation products without complication:

      - Would just have to add forward and reverse primer to the same reaction.  

      - The primers shouldn't interfere with each other.

      - I could install sequencing adaptors in the same reaction.

   - Couldn't use degenerate nucleotides (see above).

   - Would need to amplify via plasmid (see above).

Multiple copies of the target
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Extra copies of the target site give more chances for binding, which might 
improve signal.  This is something I should experiment with.

If there are multiple copies of the target, it might be hard to interpret data 
where one or more of those targets have synthesis errors.  Which target was 
actually being bound?  Motif-finding algorithms might just account for this, 
though.

Barcoding
~~~~~~~~~
Note that I do not want to barcode the DNA target:
   
- The target sequences are already quite short.

- I really want to minimize the amount of non-target DNA in the targets, 
  because I'll have to account for the fact that any DNA in the assay could be 
  bound by my proteins.

- [Hawkins2018]_ observes that synthesis errors are 10-100x more likely that 
  sequencing errors.

- But I don't want to "correct" synthesis errors.  I want to know what sequence 
  was actually being bound.

   - Note that in the topology described above with a single palindromic primer 
     between two target sites, it's most likely that the target that would be 
     amplified in not the one that would be bound.  So maybe a barcode would 
     make more sense in that situation.

- Also, I want to pick sequences based on potential for binding, not 
  information content/error correction.

- Paired-end reads can identify transcripts with sequencing errors; I might 
  just have to throw those out.

Prepare the puromycin linker
----------------------------
I'm aware of two strategies for preparing the linker.  I'm inclined to try the 
psoralen method first, because it seems easier prima facie.  But the forked 
method might be better long term.

.. note::

   All phosphoramidite prices are for 100 µmol.  Phosphoramidites are also 
   available in different kinds of bottles, for compatibility with different 
   synthesizers, but below prices don't account for this.

Psoralen linker [Kurz2001]_
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This method requires two custom synthesized oligos: one with a 3' puromycin 
modification, the other with a 5' psoralen (UV crosslinker) modification.  Both 
oligos can be synthesized by a synthesizer without manual intervention.  The 
puromycin oligo is ligated to the mRNA, then the psoralen oligo is annealed and 
crosslinked with UV light.  The psoralen oligo also primes RT.

UV crosslinking:

- Seems to be efficient.
   
   - [Naimudden2016]_ called it efficient, and they were advocating for the 
     forked linker approach.  So I'm inclined to believe them.

   - Note that the oligos are designed such that the psoralen is directly 
     across from a TpA dinucleotide, for optimal crosslinking efficiency.

- Might have to purify extra linker away.
   
   - I say "might" because I don't see what harm extra psoralen linker could 
     do.  But it's probably best to have pure reagents.

- UV light can create errors in the mRNA sequence.  However, I think UV damage 
  may be less of a concern for this application:
  
   - Crosslinking is done after translation.  At that point, only the barcode 
     sequence still matters, and it has error correction.
     
   - I'm not doing multiple rounds of mutagenesis and selection, so errors 
     won't compound. 
     
   - It's possible the proteins themselves would be damaged by UV radiation, 
     but I think this would probably be a minor concern.

Oligo sequences:

- Puromycin oligo:

   - ``5'-AAAAAAAAAACGGCTATATAAAAAAAACC-Puro-3'``,
   - 5' phosphorylated
   - Puro: Puromycin-CPG (`Glen Research`_ 20-4040, $2000)

- Psoralen oligo:

   - ``5'-Psor-TAGCCGTTTTTTTTTTAGCGGATGC-3'``
   - Psor: Psoralen C2 Phosphoramidite (`Glen Research`_ 10-1982, $195)

- Linker oligo:

   - ``5'-TTTTTTTTTTAGCGGATGC-3'``
   - Used to hold the mRNA and the puromycin oligo together during ligation.

Forked linker [Naimudden2016]_
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The forked linker is ligated to the 3' end of the mRNA.  One end of the fork 
anneals with the mRNA and acts as a primer for RT.  The other end of the fork 
contains puromycin on a flexible linker.  The basic strategy for preparing the 
forked linker is to synthesize one half of the fork, cap the end, deprotect the 
other half of the fork, the continue synthesis of this half.  Basically, it 
adds a manual deprotection step to the middle on the otherwise automated 
synthesis.  The benefit is that once this has been done, and the product 
purified, the rest of the protocol is simpler and more streamlined.

This approach has also gotten more attention recently, with method development 
papers as recent as 2016.  This makes it more likely that I could get help, if 
I needed it.

- Using a DNA synthesizer (ABI394):

   - Synthesize from 5' CCC... to 3' puromycin-CpG.
   - Deprotect branched phosphoramidite with 500 mM hydrazine hydrate.
   - Wash with pyridine/acetic acid (1:1).
   - Wash with acetonitrile.
   - Synthesize 5'-CCTG-3' from the branched phosphoramidite.
   - Elute from column using K₂CO₃.
   - Deprotect 5' acetyl group with 25% ammonium hydroxide.
   - Purify by reverse-phase HPLC.
   - Confirm by TOF-MS and gel electrophoresis.

- Product numbers (prices are for 100 µmol):

   - 5-Me-dC Brancher Phosphoramidite (`Glen Research`_ 10-1018, $205) 

     - Full name: 5'-dimethoxytrityl-N4-(O-levulinyl-6-oxyhexyl)-5-methyl-2'-deoxycytidine

   - Spacer Phosphoramidite 18 (`Glen Research`_ 10-1918, $95)

      - Full name: 18-O-Dimethoxytritylhexaethyleneglycol,1-[(2-cyanoethyl)-(N,N-diisopropyl)]-phosphoramidite

   - Fluorescein-dT (`Glen Research`_ 10-1056, $325)

      - I don't know that I actually need this, either.

   - Puromycin-CPG (`Glen Research`_ 20-4040, $2000)

In vitro transcription (IVT)
----------------------------
George recommended that I talk to Daniel Wiegand about IVT.  Daniel made the 
following recommendations:

- Use NEB PURExpress for my initial experiments.  It might be too expensive to 
  use for really large-scale experiments, but it works well.  Since my goal now 
  is just to get the assay working, it's worth using the easier and more robust 
  commercial option. 

- I should also test the in-house extracts.  If they work (e.g. not too much 
  nuclease activity, see below), they would be much cheaper.  This would be an 
  important factor for large-scale screens.

- The biggest problem I should be concerned about is that linear DNA templates 
  can get chewed up by cell extracts (e.g. various RNases, DNases).  This 
  should be less of a problem with PURExpress, but potentially more of a 
  problem with the in house extract.  I could possibly mitigate this just by 
  using more template.

- The PURExpress instructions call for 25 µL reactions, but this will produce a 
  lot of protein.  I can probably scale-down or dilute.

- The components of the PURExpress kit are His-tagged, to enable reverse 
  purification.  However, that means I *can't* use His-tags to purify my 
  proteins.

Ligation
--------
- Ligases [Bauer2017]_:

   - T4 works the best, but maybe I don't want the best.

   - T3 and human ligase 3 also work on blunt ends.
   
- Molecular crowders like PEG-6000 dramatically improve ligation efficiency 
  [Bauer2017]_.

   - But I may not want particularly efficient ligation.  I should try both.

- Sticky ends:

   - Blunt end: Presumably the best option for reducing off-target ligation.

   - Sticky end: Necessary for controlling the orientation in which ligation 
     occurs.  Depending on how the targets are designed, this control is either 
     beneficial, necessary, or unnecessary (see above).

- Phosphorylation

   - Make sure targets are dephosphorylated (so they can't polymerize).

   - Make sure cDNAs are phosphorylated.

- Quench the ligation

   - [Bauer2017]_ used proteinase K, incubated at 37°C for 30 min.  Also clears 
     the ligase off the DNA, which maybe helpful for PCR.

   - Maybe I could also just start PCR.  That would probably kill the ligases.  
     But if I'm going to be doing qPCR, I'd rather not depend on starting a 
     bunch of reactions in a consistent amount of time.


.. _Glen Research: https://www.glenresearch.com/

Methods
=======

Reverse transcribe Zif268
-------------------------
I found the wildtype Zif268 sequence on pB1H1-Zif268 from Addgene.  However, 
that sequence contains restriction sites that I might want to use, so I reverse 
translated the protein sequence myself using E. coli codons and avoiding those 
restrictions sites::

   $ sequences/reverse_translate_zif268.sh sequences/zif268/zif268.fa

Control targets
---------------
The cognate sequence for Zif268 is ``GCGTGGGCG``

As a negative control, I decided to use ``GCGAAAGCG``.  This is the negative 
control used by [Bulyk2001]_.

Make a palindromic PCR primer
-----------------------------
I wrote a script to find the skpp primer that was the most palindromic to begin 
with, had the GC content nearest to 50%, and had the longest 3' GC clamp::

   $ ./sequences/cdna_display/pick_palindromic_skpp.py

This produces a list of primers.  I chose to use::

   $ cat sequences/cdna_display/palindromic_skpp.fa
   >skpp-305-F-5
   GTATCCGAAGCTTCGGATAC

In retrospect, I decided that a palindromic primer would cause too many 
problems with PCR, and decided against using it.

Primer design
-------------
I'm going to clone into a pUC19 vector.  My plan is to do Golden Gate cloning 
with BsmBI and SapI.  This will remove most of the golden gate sites in the 
backbone, which may be useful in the future (one BsaI site will remain).  It 
will also remove 868 bp of nonessential sequence (about 30% of the plasmid), 
which should give me yields that much higher.

For the perspective of my inserts, the sticky ends will be:

- 5': ``5'-CGCG-3'`` (BsmBI)
- 3': ``3'-CGA-5'`` (SapI)

Order the Naimudden linker
--------------------------
I was able to order the DNA linker described in [Naimudden2016]_ from 
MidlandCRC.  See attached quote:

:download:`quotes/midland_20190416_63451.pdf`


