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


Methods
=======

Amplify f3-f10 --- 2020/01/21
-----------------------------
.. protocol:: 20200121_pcr.txt 20200121_dna_bead.txt 20200121_dilute_amp.txt

   See binder for corrections, e.g. use o3 instead of o13.

- It was important to do 50 µL PCR, I wouldn't have gotten enough yield with 10 
  µL.

Results
=======

Express f3 and f8 --- 2020/01/22
--------------------------------
.. protocol:: 20200122_s30_page.txt 20200122_mag_strep_page.txt

.. figure:: 20200122_s30_extract_f3_f8.svg

- I don't know why there are several green bands that are present in both 
  reactions.  These may be vaguely fluorescent proteins that are just present 
  in the S30 extract.  Next time I should include a no-template control to 
  confirm that, though.

- I don't know why I don't see any template in the +S30 lanes.  Perhaps this is 
  an indication that my template is being degraded?  The −STOP lane has a 
  yellow band stuck in the well; that could be mWasabi-repA bound to DNA, 
  although if so it doesn't bode well for the idea that the S30 extract will 
  help mWasabi-repA be well-behaved.

- I can see mWasabi expression in the +STOP reaction near the bottom of the 
  gel.  The only unique green band in the −STOP reaction is stuck in the well, 
  as discussed above.  In both cases, the level of protein expression seems 
  quite low.  

.. figure:: 20200124_s30_extract_f3_f8_streptactin_coomassie.svg

- I don't see evidence of either protein being expressed, let alone purified.  
  The gel is low quality, but I think expression is the problem.

- Don't know why the gel is so smudgy...  It looks overloaded, too.


.. todo::

   Repeat f3/f8 expression in S30 lysate, run an SDS gel, and directly image 
   mWasabi to see if the expected protein is being expressed.  Include a 
   no-template control, and maybe the provided luciferase control.

   I should also think about way to get more template.  Promega calls for 4 µg 
   per 50 µL reaction, which is a final concentration of 80 ng/µL.  I used 0.8 
   µL of ≈200 ng/µL (75 nM) template in 10 µL reactions, which is a final 
   concentration of 16 ng/µL, 20% of the recommended amount.

   The best way to get more DNA would be to do a restriction digest of plasmid.  
   Unfortunately I didn't put convenient restriction sites in these plasmids, 
   so I'd probably have to buy some enzyme that cuts the backbone in order to 
   do this.  It might just make sense to scale up the PCR.


