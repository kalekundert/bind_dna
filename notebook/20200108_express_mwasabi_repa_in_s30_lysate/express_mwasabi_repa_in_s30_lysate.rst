***********************************
Express mWasabi-repA in S-30 lysate
***********************************

As described in :expt:`20190723_confirm_cis_display_with_fluorescent_protein`, 
I think that rho factor and native RNAP may be necessary for CIS-display.  To 
test this, and to try adhering more closely to the protocol described by 
[Odegrip2004]_, I'm going to use E. coli S-30 lysate to express mWasabi-repA 
fusions.

Considerations
==============

Meaning of "S-30"
-----------------
S-30 refers to the "supernatant of the 30,000g centrifugation" (or maybe 
"sedimentation at 30,000g", but you get the point) after lysing the cells by 
French press [Lesley1995]_.  It is a method for preparing cell lysates, and can 
be used for many different strains.

Genotype requirements
---------------------
According to [Lesley1995]_ (Table 1, pg 268), the following genotypes are 
required in order to express genes from linear templates:

- recBC or recD: These genes encode 3 proteins which assemble into a helicase.  
  The helicase is active against linear DNA, and probably interferes with 
  transcription from linear templates.
- sbc: Exonuclease that is somehow associated with the recBCD complex 
  `[EcoliWiki] 
  <https://ecoliwiki.org/colipedia/index.php/sbcC:Gene_Product(s)#cite_note-LIB:EcoGene-4>`_.

In addition, the following genotypes are recommended to improve protein yield:

- hsdS: Reduce restriction digestion.
- ompT: Cleaves at bibasic residues [Lesley1995]_.
- lon: Heat shock protease [Lesley1995]_.

The SL119 strain is recommended for linear templates.  The A19 strain is not 
mentioned.  As discussed above, I think A19 is a relatively wildtype strain and 
probably has all of the above exonucleases and proteases.

A19 genotype
------------
The A19 strain is a common choice for in vitro protein expression.  Dan Weigand 
has a stock of A19, and offered to show me how to prepare an extract from it.  
The `A19 genotype <https://cgsc2.biology.yale.edu/Strain.php?ID=7376>`_ is 
described on the CGSC website.  The strain was discovered during a screen of 
chemically mutagenized E. coli for lack of RNase activity [Gesteland1966]_.  
The particular mutation causing this lack of activity is termed "rna-19", and 
is located in the gene for RNase A (note that older papers refer to RNase A as 
RNase I).  Other mutations in the A19 strain seem to be inherited from the 
parent strain (AB301, which I believe to be related to DH10B, and therefore in 
the K-12 family) and do not seem particularly relevant to protein expression.

Presumably A19 became a popular strain for in vitro transcription/translation 
due to its lack of RNase activity.  Other than that, it seems pretty wildtype.

SL119 genotype
--------------
SL119 is the strain used by Promega to create their "E. coli S30 Extract System 
for Linear Templates" (Promega L1030) kit.  It is also the strain recommended 
for linear templates by [Lesley1995]_.  It is derived from the B strain (e.g.  
in the same family as BL21, not in the same family as A19) and has the 
following genotype [Lesley1995]_:

- hsdS: Restriction enzymes
- gal: Lower background for galactokinase assays (I think).
- ompT: Protease
- recD:Tn10: Helicase; the \*:Tn10 notation indicates that the gene, recD in 
  this case, is disabled due to insertion of the Tn10 transposon.

Overall it seems like I really want SL119 and not A19.  I'll ask Dan if he has 
that strain, otherwise I'll just order the extract from Promega.

.. update:: 2020/01/21

   I ordered from Promega.

DNA purification
----------------
Promega makes the following recommendation for PCR-amplified templates:

   Avoid contaminating the S30 reaction with the wrong PCR product or primer 
   dimers. If agarose gel analysis indicates that your PCR produced a unique 
   band, primer dimers can be removed by ethanol precipitation with sodium 
   acetate.  Otherwise, PCR-amplified DNA should be gel-purified before using 
   it in the E. coli S30 Extract System for Linear Templates.

My amplicons looks pretty clean (with o3 and o68).  I think a bead purification 
of silica column would also serve to get rid of excess primer, so I'm going to 
try those before doing an ethanol precipitation.




