************************
Linker with click handle
************************

In order to do cDNA display, I'll need to attach a reverse transcription (RT) 
primer in the linker used by [Doshi2014]_.  The goal of this experiment is to 
find a way of doing this that doesn't interfere with the display itself.

Considerations
==============

Attachment chemistry
--------------------
I found previously (:expt:`57`) the azide/DBCO click chemistry is an easy and 
robust method for attaching oligos together.

- Details on the azide modification:

    At IDT, we incorporate azides into oligos after synthesis, using NHS ester 
    chemistry.  An amino modification, added during oligo synthesis, is used in 
    the NHS ester reaction. For example, a 5' azide is attached to the oligo 
    sequence through a `5' Amino Modifier C6 <5AmMC6>`_, while a 3' azide is 
    attached through a `3' Amino Modifier <3AmMO>`_. An internal azide is 
    attached through the `Amino Modifier C6 dT <iAmMC6T>`_, which results in 
    the incorporation of a dT base at that location in the oligo sequence.

- Will two click-attached nucleotides be able to base pair?

  - Looking at the crystal structure of duplex DNA and the chemical structure 
    of the internal C6 dT nucleotide, it's very unlikely that the azide 
    modification will interfere with base pairing.

  - C-C bond counts for between relevant atoms in various modified nucleotides:

    - alkyne to 5' P (5DBCON_): 21
    - alkyne to amide N (5DBCON_): 10
    - azide to amide N (5AzideN_): 7
    - azide to alkyne (triazole): 1-2
    - dT methyl to free amine (iAmMC6T_): 10
      - Used for 5DBCON_ and 5AzideN_.

  - Total C-C bond counts for different click pairs:

    - 5AzideN_ (dT methyl) to 5DBCON_ (dT methyl): 10+7+1+10+10 = 38
    - 5AzideN_ (dT methyl) to 5DBCOTEG_ (5' P): 10+7+1+21 = 39

  - Distances in B-form DNA (measured in PDB ID: 1BNA)

    - Duplex: /iAzideN/A + /5DBCON/A

      - dT methyl to dT methyl: ≈8.4Å
      - Modifications offset by 1
      - This sequence doesn't actually appear in 1BNA, so I just got as close 
        as I could.

    - Duplex: AN/iAzideN/ + /5DBCON/NA

      - dT methyl to dT methyl: 11.2Å
      - Modifications offset by 2

    - Duplex: /iAzideN/ + /5DBCOTEG/A

      - dT methyl to 5' P: 14.1Å

    - Duplex /iAzideN/C + /5DBCOTEG/GA

      - dT methyl to 5' P: 12.4Å
      - Modifications offset by 1 nt.
      - This occurs with the Nt.BstNBI site

    - Duplex /iAzideN/NN + /5DBCOTEG/NNA

      - dT methyl to 5' P: 12.5Å
      - Modifications offset by 2 nt.

    - Duplex /iAzideN/NNN + /5DBCOTEG/NNNA

      - dT methyl to 5' P: 13.3Å
      - Modifications offset by 3 nt.

  - It doesn't really seem to matter if I use 5DBCON_ or 5DBCOTEG_; both are 
    about the same length and are easily long enough to bridge small offsets.  
    I can instead decide which to use based on (i) price and (ii) if there's a 
    dA for 5DBCON_ to base pair with.

  - The 1 nt offset that occurs in the Nt.BstNBI site actually seems about 
    optimal for allowing the attached oligos to base pair.  The distance is 
    minimized, and the straight-line path does not go through any other 
    atoms.

  - It seems very likely that two click-modified nucleotides will be able to 
    base pair.

Click timing
------------
- Before ligation:

  - This is ideal because I can use very high concentrations of oligo, and I 
    can be sure (by doing gel purification, or just observing the reaction go 
    to completion) that every linker has all three arms.

  - Greater chance that the attached oligo interferes with a downstream step.

- During/after ligation:

  - Minimize chance that ligation fails.

  - In the case of LR-P, can purify only mRNAs with puromycin.

- After expression:

  - Only an option for LP-R.

  - Very unlikely that oligo interferes with ligation or display.

  - Relatively high likelihood that unattached oligos could prime RT.

  - Proteins might denature during long click reaction.

T4 RNA ligase 1
---------------
I was reading about T4 RNA ligase 1 [England1977]_ [Wang2006]_ to see if it 
would be affected by the presence of an click-attached oligo.  Below are some 
things I learned:

- The basic mechanism is: [Wang2006]_

  - RNA ligase reacts with ATP to form a covalent ligase-(lysyl-N)–AMP 
    intermediate plus pyrophosphate.
  - AMP is transferred from ligase-adenylate to the 5′-PO4 RNA end to form an 
    RNA–adenylate intermediate.
  - ligase catalyzes attack by an RNA 3′-OH on the RNA–adenylate to seal the 
    two ends via a phosphodiester bond and release AMP.

- Instead of ordering phosphate-modified oligos, I could order 
  5'-adenosine-pyrophosphate-modified oligos.  These would react even in the 
  absence of ATP.  Not sure if there's any practical advantage to this, though, 
  and this modification is very expensive.  For more info, see `5' Adenylation 
  <5rApp>`_.

- Crystal structure: 2C5U [ElOmari2006]_

Puromycin arm length
--------------------
Translation will stop at the RNA/DNA junction.  For there, according to 
[Keefe2001], the linker should be very nearly 30 nt (including the puromycin, 
which mimics an A).  [Keefe2001] claims that large deviations from 30 nt will 
cause dramatic decreases in yield.  

- The [Doshi2014]_ linker (splint-ligation)::

    /5Phos/GCAAAAAAAAAAAAAAAAAAA/iSp18//iSp18/ACC/3Puro/
    length=27

- The [Reyes2021]_ linker (Y-ligation)::

    /5Phos/CCCTTCACCTGATCCGCTGAAAAAAAAAAAAAAAAAA/iSp18//iSp18//iFluorT//iSp18/CC/3Puro/
    length=44

- The [Barendt2013]_ linker (splint-ligation)::

    AAAAAAAAAAAAAAAAAAAAAAAAAAACC/3Puro/
    length=30

- The [Keefe2001]_ linkers (splint ligation)::

    AAAAAAAAAAAAAAAAAAAAAAAAAAACC/3Puro/
    length=30

    AAAAAAAAAAAAAAAAAAAAA/iSp9//iSp9//iSp9/ACC/3Puro/
    length=28

Nicking endonuclease
--------------------
It might be useful to include a nicking endonuclease site in the RT primer, so 
that I can protect its 3' end from ligation and later free it for 
transcription.  Here are some data on all the nicking endonucleases that NEB 
sells:

.. datatable:: nicking_endonucleases.xlsx

- Nt.BstNBI and Nt.BsmAI seems like the two best candidates.

  - The advantage of Nt.BstNBI is that it is a natural enzyme, which means that 
    it should have absolutely no double-stranded cutting activity.  The 
    engineered (i.e. non-natural) enzymes still have the active site to cut the 
    opposite strand, it's just been inactivated with mutations.

  - The advantage of Nt.BsmAI is that it is shorter.

  - Both contain dT, which can be used to install the azide.

Oligo topology
--------------
There are two possible ways to attach two oligos such that the resulting 
product has a free 5' end (for ligation), a 3' puromycin (for mRNA display), 
and a free 3' end (for RT):

.. figure:: linker_topologies.svg
   :align: center

====  ====================================  ==================================
Name  Oligo 1                               Oligo 2
----  ------------------------------------  ----------------------------------
LR-P  5'-[ligation]-azide-[RT]-3'           5'-DBCO-[puromycin]-3'
LP-R  5'-[ligation]-azide-[puromycin]-3'    5'-DBCO-[RT]-3'
====  ====================================  ==================================

LR-P
~~~~
.. datatable:: lr_p_oligos.xlsx

::

  LR: /5Phos/[splint'(with /iAzideN/)][RT'][hairpin][RT]

  P: /5DBCON/AAAAAAAAAAAAAAAAAAAAAAAA/iSp18//iSp18/ACC/3Puro/

- I already have several puromycin arms, e.g. o125

- Where to attach the puromycin arm:
  
  - As close as possible to the 5' end of the ligation arm:

    - Keeps the number of nucleotides between the puromycin and the ribosome 
      about the same.

    - The oligo-dT purifications should be unaffected, because there will still 
      be a ≈30 nt oligo-dA arm.
    
    - Might interfere with ligation.
      
      - Might be able to avoid this problem by using 5' adenylation instead of 
        5' phosphorylation.  T4 RNA ligase seems pretty promiscuous w.r.t. the 
        5' adenylated substrate [England1977]_ , but it's possible that the 
        adenylation step itself might be more sensitive.
        
      - This problem could also be side-stepped by doing the click reaction 
        after the ligation reaction.

    - I can use a long primer (e.g. Nt.BstNBI), since it doesn't affect the 
      length of anything else.

  - ≈10-15 nt from the 5' end of the ligation arm:

    - Unlikely to interfere with ligation.

    - May still interfere with reverse transcription.

    - May need to shorten the oligo by ≈10-15 nt, to keep the distance between 
      the puromycin and the ribosome near what it should be.  The oligo-dT bead 
      purification might stop working, though.

  - Behind the RT primer hairpin:

    - Unlikely to interfere with ligation or reverse transcription.

    - Most likely to need length optimization.  The double-stranded segment 
      will affect flexibility, so the 30 nt recommendation by [Keefe2001]_ may 
      not apply.

    - Quite similar to LP-R.

- Splint complementarity:

  - The [Doshi2014]_ linker is mostly dA, so the splint is mostly dT.

  - I need oligo-dA for bead purifications, but that will be in the P-arm, so I 
    can chose any sequence I want for the splint (as long as it includes a dT 
    to attach the azide to).

  - I probably don't want to use oligo-dA, because I'd rather the beads bind to 
    the puromycin arm.  This would make it possible to do the ligation and 
    click reactions simultaneously.

  - Length: [Keefe2001]_ recommends 10 nt, [Doshi2014]_ uses 15 nt.
    
  - If attaching the P-arm before the hairpin:

    - In this case, the length of the splint doesn't matter, so I can just use 
      a full SR primer.

    - I'll put the reverse complement of the primer in the oligo, so that the 
      primer itself can be used to amplify the gene.

    - I want the splint sequence in the LR-oligo (i.e. the reverse complement 
      of the SR primer) to have lots of dT nucleotides, so that I can vary 
      where the P-arm is attached without having to reorder the splint.  So I 
      want an SR primer with lots of dA nucleotides.

      - SR067

      - I sorted all the SR primers by these metrics:

        - ``a_islands``: ``re.sub('A+', ':', x).count(':')``
        - ``a_count``: ``df.sequence.str.count('A')``

      - Of the primers with the most A-islands, SR067 and SR025 are the only 
        ones with no two adjacent A's.  I like that, because it means the A's 
        are well spaced.  I prefer SR067 because it has more A's near the 5' 
        end.

      - Would be convenient if there were a G in the last 5-7 nt, to pair with 
        /3ddC/.  That doesn't really matter, though.

    - This probably isn't really relevant for ssDNA, but there are 10.5 bp per 
      helical turn in dsDNA.  So I might look for dT nucleotides spaced by 5 to 
      see if it matters which side of the helix the P-arm is attached to.

  - If attaching the P-arm after the hairpin:

    - I should probably go back to the oligo-dA splint complementarity 
      sequence, because the P-arm probably can't be long enough to bind the 
      oligo-dT beads itself.

    - That said, I might as well still try ordering 1 linker with more-or-less 
      the same sequence as the before-the-hairpin designs, just to see.

- Need a hairpin:

  - [Kannan2007]_: ``gcGCAgc``

  - [Vallone1999]_:

    - Experimentally characterize 28 16 nt hairpins, all with 6-bp duplex, 4-nt 
      turn, and identical first 5 bp.

    - None of these hairpins are belong to the families identified in 
      [Nakano2002]_.  Frankly, they all seem like pretty awful hairpins.

    - Melting temperatures from 32.4-60.5°C.

      - Could be useful to titrate melting temperature.
      - NEB calls for the MMLV RT reaction to occur at 37-42°C.

  - [Nakano2002]_:

    - Screen for stable tetraloops

    - Use 12 bp stem in screen, designed to avoid slipped pairings and to 
      include some cloning sequences.

    - Four families identified:

      - ``cGNNAg``
      - ``cGNABg``
      - ``cCNNGg``
      - ``gCNNGc``

    - "For d(cGNNAg) and d(cGNABg), a CG closing base pair was strongly 
      preferred over a GC, with ΔΔG°37 ≈ 2 kcal/mol."

  - [VanDongen1997]_:

    - NMR structure of an ``aGTTAt`` tetraloop.

    - Very likely that either dT could be azide-modified without clashing with 
      anything.  I'd probably choose to modify the first, since it's not even 
      stacking with anything.

  - Seems like there are plenty of options, and this probably isn't a very 
    important parameter.

  - I'm going to use ``cGTTAg`` to start with:

    - Known sequence [Nakano2002]_.
    - Known structure [VanDongen1997]_.
    - Can attach modifications (e.g. azides) to dT nucleotides.
    - Nt.BstNBI has CG base pair right before the hairpin, although that won't 
      matter if I decide to extend the duplex a bit.

- How long should the primer be?

  - Probably doesn't need to be very long, by virtue of being unimolecular.

  - I should test this empirically:
      
    - Buy cheap primers with different hairpins.
    - Kinase, splint-ligase.
    - MMLV RT
    - RNase
    - Gel

  - Random hexamers are used for RT, so it's unlikely that the RT needs more 
    duplex than that.

  - I'm going to use Nt.BstNBI for the designs where the P-arm is attached 
    before the duplex.

    - The length of the primer is not important for these designs.

    - This restriction site allows a 3' cap to be removed, which may be 
      necessary (see discussion below).

    - It might be necessary to extend the duplex beyond just the recognition 
      sequence, because not all restriction enzymes cleave well at DNA ends.  
      But I'm not going to worry about that for now.  The goal of this 
      experiment is to find linkers that are compatible with ligation and 
      puromycin coupling.  If those steps work, I'll be able to optimize the 
      length of the duplex later.

- Capping the 3' end:

  - I'm not sure if this is necessary.

  - IDT requires the 3' end to be blocked (with either /3SpC3/ or /3ddC/) when 
    ordering an oligo with the 5rApp_ modification, to avoid circularization.  
    But I don't know if this concern is specific to 5rApp_, or if it's 
    something that needs to be considered whenever using T4 RNA ligase.  The 
    former might be the case if 5rApp_ is a good enough leaving group for the 
    oligo to circularize on it's own, without an enzyme.  The latter might be 
    true if it's T4 RNA ligase that's catalyzing the circularization in both 
    cases.

  - The linker in [Doshi2014]_ naturally has it's 3' end protected by the 
    puromycin modification.

  - The safest thing to do for now is to include the cap.  I can experiment 
    with removing it once I'm trying to get the RT step working.

- Puromycin arm:

  - 5DBCON_: has a normal 5' end (the DBCO is attached to the methyl group in 
    the thymidine base), so there will be a free 5'-OH.  However, this end 
    shouldn't be a substrate for lambda exonuclease or T4 RNA ligase, since it 
    isn't phosphorylated.

  - 5DBCOTEG_: the DBCO is basically attached to iSp9_, and the DBCO moiety 
    itself has 9 heavy atoms leading up to the alkyne, so this effectively 
    includes an iSp18_ spacer.  The azide is also on a fairly long liker.

- Ordering:

  - I don't really need to worry about the hairpin/primer yet, because those 
    decisions are unlikely to affect the puromycin coupling reaction.

  - So I should just order something reasonable.  At the same time, I can order 
    a bunch of cheaper oligos to see which hairpin/primer sequences work best 
    for RT.

LP-R
~~~~
.. datatable:: lp_r_oligos.xlsx

- This topology will stay much more similar to the [Doshi2014]_ linker.

- To avoid interfering with the splint, the double-stranded region has to end 
  >10 nt from the 5' end of the ligation arm.  Of course, this won't matter 
  if the splint is attached after ligation.

  - The [Doshi2014]_ splint is 15 nt, but 10 nt splints have been used in other 
    papers.

  - Without increasing the length of the LP oligo, the RT oligo will have to be 
    quite short (or added after ligation).

    - If the RT oligo includes a nicking restriction site, there won't really 
      be room to extend the duplex any more than for the site itself.  I don't 
      know if the restriction enzyme will work like that.

  - There are about 35 C-C bonds between the dT methyl (attached to the azide) 
    and the 5' phosphate (attached to the DBCO-TEG).  The distance between 
    these two atoms in B-DNA (PDB: 1BNA) is 14.2Å.  Therefore, it seems likely 
    that the 5' nucleotide in the RT arm can base-pair with the /iAzideN/ 
    nucleotide.

- The splint is phosphorylated.  That means I need to make it doesn't get 
  ligated to the primer.

  - Attach the primer after ligation.
  - Block the 3' end of the primer.

- The increased rigidity of the double-stranded region might interfere with 
  the ability of the puromycin to enter the ribosome, but I think this is 
  unlikely to be a problem as long as the spacers come afterward.

- Two basic designs:

  - Just add the azide, changing nothing else.

  - Add the azide with a nicking endonuclease site.  Extend the length of the 
    linker slightly.  This will allow me to order RT arms with blocked 3' ends, 
    allowing the RT arm to be attached before ligation.



.. _5AmMC6: https://www.idtdna.com/site/Catalog/Modifications/Product/1082
.. _3AmMO: https://www.idtdna.com/Site/Catalog/Modifications/Product/3299
.. _iAmMC6T: https://www.idtdna.com/Site/Catalog/Modifications/Product/1388
.. _5DBCON: https://www.idtdna.com/pages/education/decoded/article/need-a-non-standard-modification
.. _5DBCOTEG: https://www.idtdna.com/pages/education/decoded/article/need-a-non-standard-modification
.. _5AzideN: https://www.idtdna.com/Site/Catalog/Modifications/Product/3770
.. _5rApp: https://www.idtdna.com/Site/Catalog/Modifications/Product/2342
.. _iSp9: https://www.idtdna.com/site/Catalog/Modifications/Product/1391
.. _iSp18: https://www.idtdna.com/site/Catalog/Modifications/Product/1393
