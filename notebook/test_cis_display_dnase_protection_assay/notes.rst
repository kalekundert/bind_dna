****************************************
Test CIS-display, DNase-protection assay
****************************************

I have confirmed that Zif268 binds its TGG target (:expt:`35`) and shown that 
CIS-display may be successfully linking Zif268 to the DNA encoding it (the data 
can be explained in other ways, too, but is not inconsistent with CIS-display 
working as intended; :expt:`36`).  The next step is to test whether or not this 
CIS-display construct can be used to detect DNA binding between Zif268 and its 
target.

I want to begin with the DNase protection assay suggested by Kettner, because I 
like that its 1-body nature eliminates concerns about the kinetics of proteins 
finding their targets.  But I hope to pursue the ligation based assay in 
parallel.

The two most important considerations for this assay are:

- Does binding provide adequate protection from DNase treatment?
- Do *cis* interactions outweigh *trans* interactions?

My plan for testing both of these considerations involves two protein/DNA 
pairs: ``Aa`` and ``Bb``, where I'm referring to the proteins as ``A`` and 
``B`` and to the corresponding DNA targets are ``a`` and ``b``

- Using qPCR, I first need to show that ``Aa`` and ``Bb`` produce high signal, 
  while ``Ab`` and ``Ba`` produce low signal.  This would indicate that the 
  ``A``/``B`` binding does in fact provide protection from DNase treatment.
 
- Next, I can test that the 1:1 mixture ``Ab + Ba`` does not produce signal.  
  The proteins in this mixture could only provide DNase protection via trans 
  interactions, so the lack of signal indicates a condition where trans binding 
  does not occur.

Note that it isn't strictly necessary for ``Bb`` to be a functional pair.  I 
could test the first point using just Zif268 (``A``) and the TGG (``a``) and 
AAA (``b``) target/non-target sequences.  I could also test the second point 
using Zif268 (``A``) and repA-without-Zif268 (``B``).

Considerations
==============

Protein/DNA target pairs
------------------------
As outlined above, I need orthogonal protein/DNA target pairs to test the DNase 
protection assay.  Obviously Zif268/TGG will be one pair, but there are 
multiple possibilities for the second, as discussed below.

For now, I've decided to use REDV and LRHN from [Bulyk2001]_ and the ΔZif268 
control.  I think this will give me a nice set of options for getting the assay 
working and confirming that the dynamic range is compatible with testing my 
designs.

Nonfunctional
~~~~~~~~~~~~~
I really only need one functional protein/DNA target pair to be able to test 
that the assay works and to distinguish between cis and trans interactions.  In 
fact, the data might be easier to interpret if there's only one positive 
interaction that can occur.  This might be the way to go to get the assay set 
up, but I'll also need controls that better mimic the dynamic range I can 
expect when testing designs.

- GFP or any other non-Zn-finger

   - Obviously these proteins would have no activity.

   - The folding and display properties of these proteins may differ

- Just repA, i.e. ΔZif268.

   - Seems like this would be a pretty intuitive negative control.  Elegantly 
     controls for any repA binding, too.  (I mean, so would the other ideas, 
     but this would really let me know if that was happening, because nothing 
     else is going on.)

- AAAA

   - I could always make my own negative control

   - I wouldn't know a priori that my mutant folds correctly and doesn't have 
     any residual activity.

   - See [ElrodErickson1998]_ for affinity data for 5 alanine mutants in the 
     1st finger.  The strongest effects are for R18A (-1) and R24A (6), and are 
     on the order of 100x.

Orthogonal---Artificial
~~~~~~~~~~~~~~~~~~~~~~~
- RGPD/GCG, REDV/GCG, LRHN/TAT, KASN/NNT [Bulyk2001]_

   - These mutants have the best characterized specifities w.r.t. wildtype of 
     any that I've found.

   - The positions these sequences correspond to is given in Figure 1D.

   - The mutants came from a screen that was done specifically to test this 
     method.  However, most (if not all) of the sequences seems to have been 
     found previously by [Choo1994]_.  [Choo1994]_ also describes a number of 
     other sequences that were selected in the same manner and may be good 
     controls, but have not been characterized as extensively as the sequences 
     in this paper.

   - The screen does not appear to include any counter-selection that would 
     encourage specificity, e.g. against the wildtype target sequence, although 
     the specificities are measured.

   - These mutants may not be the best in terms of absolute orthogonality, 
     since they do not differ much from wildtype Zif268.  But they would tell 
     me realistically if my dynamic range is appropriate, because the designs I 
     want to test would probably have similar numbers of mutations (focusing on 
     a similar part of the binding site).

   - "The mutant RGPD was selected from the phage display library by using the 
     DNA sequences GCG and GCT, whereas the mutant REDV was selected by using 
     GCG and GTG. The microarray-binding results not only verify binding to 
     these respective sequences but also indicate that RGPD binds fairly well 
     to CCG and to GCT (although RGPD was not recovered from a selection by 
     using CCG)."

     "The mutants LRHN and KASN had been isolated repeatedly after independent 
     sets of in vitro selections by using many different 3-bp binding sites for 
     the second finger (ACT, AAA, TTT, CCT, CTT, TTC, AGT, CGA, CAT, AGA, AGC, 
     and AAT).  [The authors therefore considered LRHN and KASN to be examples 
     of zinc fingers with poorly characterized sequence specificity. -KBK]"
     
     "Although LRHN was isolated in all of these selections, microarray-binding 
     experiments revealed that this variant is highly specific for the DNA 
     sequence containing the 3-bp binding site TAT (Table 2C of supplemental 
     data). Moreover, the LRHN–TAT complex is almost as tight as the wild-type 
     Zif268–DNA complex in these experiments."
     
     "Meanwhile, microarray-binding experiments showed KASN to be fairly 
     nonspecific, with binding to a number of DNA sequences consistent with a 
     central 3-bp consensus (A/C/T)NT (Table 2D of supplemental data)."

   - REDV seems like the most specific of these sequences.  LRHN may be a good 
     candidate, too.

   - REDV was screened to bind GCG/GTG; and has affinity for both.  Even 
     knowing that REDV is orthogonal to Zif268, I'd rather a finger that was 
     optimized to bind a single sequence.

- TATA box, p53 site, nuclear receptor element (NRE) site [Greisman1997]_

   - Nothing really special to recommend, just examples of Zn-fingers that have 
     been selected to bind particular targets.

- Aart
   
   - Used by [Zykovich2009]_

   - 6-finger protein designed to recognize ANN triplets.

   - Would rather a Zif268 mutant, to minimize the differences between my two 
     controls, even though the assay should work with any two orthogonal DNA 
     binding proteins.

- DSNR, QGSR, RADR [Rebar1994]_, [ElrodErickson1998]_

   - These were hits from [Rebar1994]_, which I think was the first selection 
     of zinc-fingers with altered specificities.

   - I can't get access to [Rebar1994]_ (old Science paper...), but the 
     sequences are also described in [ElrodErickson1998]_, where they are 
     crystallized.

   - These designs are named following the convention that NNNN refers to 
     positions -1, 2, 3, and 6, counting from the beginning of the α-helix.  
     For example, the wildtype sequence of the first finger of Zif268 is RDER.  

   - DSNR:

      - 1st finger: bind GACC
      - A structure is also reported for DSNR binding the wildtype GCGT 
        sequence, which it "binds less tightly".  Still, there is clearly some 
        overlap, so this wouldn't be a good orthogonal control.

   - QGSR:

      - 1st finger: GCAC
      - Wildtype sequence used as competitor during screen.
      - A structure is reported for wildtype Zif268 binding the "less 
        favorable" GCAC site.  Not clear if this would cause a problem.

   - RADR:

      - 1st finger: GCAC
      - Nonspecific DNA used as competitor during screen.
      - Mutated residues form fairly nonspecific (e.g. DNA backbone or protein) 
        contacts.
      - RADR still binds "very tightly" to wildtype GCGT sequence, so it would 
        clearly not be a good orthogonal control.

- AZP4 [Sera2002]_, [Nakata2012]_

   - Artificial Zn-finger Protein (AZP)

   - These are designed without selection using code-based rules.  I don't 
     trust that process to produce specific binders.

Orthogonal---natural
~~~~~~~~~~~~~~~~~~~~
I could also use other natural zinc fingers.  The advantage of these proteins 
is that I can expect them to be well-folded and highly orthogonal to Zif268.  I 
really just need to figure out which is the most well-behaved.  The drawbacks 
are that they introduce more differences that a Zif268 mutant would, and are 
not necessarily optimized for specificity.

- GCN4

- TFIIIA

- ADR1



Methods
=======

Cloning --- 2019/10/18
----------------------
.. protocol:: 20191018_pcr.txt 20191018_golden_gate_n.txt 20191018_golden_gate_c.txt 20191018_golden_gate_c_null.txt

   Summary of all cloning reactions: :download:`sequences/reactions.xlsx`

   ***

   .. figure:: 20191018_golden_gate_pcr.svg

   All of the PCR reactions worked pretty well.  The only one with some 
   problems was L, which was used only to make 79.  We'll see if I have any 
   trouble with that reaction in particular.

   ***
   ***
   ***

   For all the the N-terminal inserts, I got a mix of big and small colonies.  
   I'm assuming that the big colonies are more likely to be correct, but I 
   don't know what the small colonies represent.

.. protocol:: 20191021_pcr.txt

   In retrospect, I don't think these were very good primers for checking the 
   N-terminal inserts.  Those inserts only add 0.4 kb to the amplicon.  It 
   probably would've been smarter to use a reverse primer in repA or something.

   ***

   .. figure:: 20191021_colony_pcr.svg

   Note that I accidentally switched 63 and 69.  The bands for 77-81 didn't 
   seem to run straight, so they're hard for me to interpret.  

Results
=======

qPCR --- 2019/11/01
-------------------
Before doing PURExpress reactions, I want to make sure that I can tell the 
difference between digested and undigested DNA with qPCR.  The basic idea is to 
incubate linearized plasmid (59) with and without exonuclease.  



After this, once I start doing qPCR, I should include free Zif268 (34) as a 
control.
