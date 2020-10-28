********************************
Confirm dCas9-PCV2 concentration
********************************

In :expt:`68`, I became worried that the dCas9 aliquots I got from Jorge might 
be ≈10x less concentrated than I thought they were.  To check this, I measured 
the concentration of the aliquot by Nanodrop.  I wasn't too concerned about 
accuracy, because I just wanted to know if I was in the right order of 
magnitude.

Considerations
==============

Buffer
------
I wasn't sure what buffer the protein was in, in order to blank the Nanodrop.  
I looked up the purification protocol I got from Jorge and assumed that the 
protein was in the storage buffer from that protocol:

- 50 mM HEPES, pH 7.5
- 100 mM NaCl
- 1.5 mM EDTA
- 20% glycerol
- 1 mM DTT

I didn't have all the chemicals to make that exact buffer, so I just did the 
best I could with the chemicals I had on hand:

- 100 mM NaCl
- 1.5 mM EDTA
- 20% glycerol

I also confirmed that HEPES is not expected to absorb at 280 nm: 
:download:`buffer_a280.pdf`

Extinction coefficient
----------------------
I used the ExPASy ProtParams tool to predict the extinction coefficient for the 
dCas9-PCV2 fusion.  I used the prediction that assumes all cysteines are 
engaged in disulfide bonds::

  Number of amino acids: 1538

  Molecular weight: 176760.70

  Theoretical pI: 9.14

  Amino acid composition: ￼CSV format
  Ala (A)  82	  5.3%
  Arg (R)  87	  5.7%
  Asn (N)  77	  5.0%
  Asp (D) 102	  6.6%
  Cys (C)   5	  0.3%
  Gln (Q)  59	  3.8%
  Glu (E) 121	  7.9%
  Gly (G)  92	  6.0%
  His (H)  40	  2.6%
  Ile (I)  97	  6.3%
  Leu (L) 155	 10.1%
  Lys (K) 167	 10.9%
  Met (M)  24	  1.6%
  Phe (F)  69	  4.5%
  Pro (P)  45	  2.9%
  Ser (S)  99	  6.4%
  Thr (T)  72	  4.7%
  Trp (W)   9	  0.6%
  Tyr (Y)  58	  3.8%
  Val (V)  78	  5.1%
  Pyl (O)   0	  0.0%
  Sec (U)   0	  0.0%

   (B)   0	  0.0%
   (Z)   0	  0.0%
   (X)   0	  0.0%

  Total number of negatively charged residues (Asp + Glu): 223
  Total number of positively charged residues (Arg + Lys): 254

  Atomic composition:

  Carbon      C	      7899
  Hydrogen    H	     12569
  Nitrogen    N	      2191
  Oxygen      O	      2350
  Sulfur      S	        29

  Formula: C7899H12569N2191O2350S29
  Total number of atoms: 25038

  Extinction coefficients:

  Extinction coefficients are in units of  M-1 cm-1, at 280 nm measured in water.

  Ext. coefficient   136170
  Abs 0.1% (=1 g/l)   0.770, assuming all pairs of Cys residues form cystines

  Ext. coefficient   135920
  Abs 0.1% (=1 g/l)   0.769, assuming all Cys residues are reduced

  Estimated half-life:

  The N-terminal of the sequence considered is M (Met).

  The estimated half-life is: 30 hours (mammalian reticulocytes, in vitro).
                              >20 hours (yeast, in vivo).
                              >10 hours (Escherichia coli, in vivo).

  Instability index:

  The instability index (II) is computed to be 41.00
  This classifies the protein as unstable.

  Aliphatic index: 83.94

  Grand average of hydropathicity (GRAVY): -0.648


Results
=======

Nanodrop --- 2020/10/16
-----------------------
I measured 2 µL of undiluted aliquot using the Nanodrop2000.  I got the 
following results:

- Concentration: 3.419 mg/mL
- A280 (adjusted to 1 cm path length): 2.634
- 260/280: 0.71

Thermo states that the A280 assay with the NanoDrop 2000 is accurate within a 
range of 0.1–400 mg/mL: :download:`protein_measurements.pdf`.  This value is 
right in the middle of that range, so I have reason to think it's accurate.

Because I was a little confused by the interface for providing the extinction 
coefficient to the NanoDrop program, I used the measured absorbance value to 
calculate the protein concentration myself as a sanity check.  The following is 
Beer's Law, which gives absorbance (:math:`A`) as a function of path length 
(:math:`b`), concentration (:math:`c`), and the molar extinction coefficient 
(:math:`\epsilon`):

.. math::

  A = b c \epsilon

Rearranging for concentration and plugging in the appropriate values:

.. math::

  c &= \frac{A}{b \epsilon} \\
    \\
    &= \frac{2.634}{\pu{1 cm} \cdot \pu{136170 M-1 cm-1}} \times
       \frac{\pu{10^6 µM}}{\pu{1 M}}\\
    \\
    &= \pu{19.34 µM} \\
    \\
    &= \frac{\pu{19.34 µmol}}{\pu{1 L}} \times
       \frac{\pu{176760.70 µg}}{\pu{1 µmol}} \times
       \frac{\pu{1 mg}}{\pu{1000 µg}} \times
       \frac{\pu{1 L}}{\pu{1000 mL}} \\
    \\
    &= \pu{3.419 mg/mL}

The dCas9-PCV2 fusion is labeled as being 22 µM, and this measurement is very 
much consistent with that.

Because there really seems to be an order of magnitude more DNA than protein in 
my reactions (on a molar basis), my next instinct is to check my math.  As 
detailed in the footnotes of my HUH-tagging protocol, Invitrogen recommends 
loading 250 ng protein per band: :download:`bolt_quick_ref.pdf`

That corresponds to:

.. math::

  \pu{250 ng} \times
    \frac{\pu{1 nmol}}{\pu{176760.7 ng}} \times
    \frac{\pu{1000 pmol}}{\pu{1 nmol}}
  = \pu{1.4 pmol}

I've been using 1 pmol per reaction, which is comparable to that.  Needless to 
say, this does not seem like enough.  Putting that aside, though, here is the 
calculation for how to add a 1:1 molar ratio of DNA:

.. math::

  \pu{1.4 pmol} \times
    \frac{\pu{260299.81 pg}}{\pu{1 pmol}} \times
    \frac{\pu{1 ng}}{\pu{1000 pg}}
  = \pu{368 ng}

My DNA has been about 150 ng/µL, and I've been adding about 2 µL per reaction, 
so that all seems right.  I don't see any mistakes.

My thought now is that if Jorge and I both used the NanoDrop to measure the 
protein concentration, and there's some reason why the NanoDrop would 
overestimate the concentration, that would explain the results I've observed.  
The way to test this hypothesis would be to use a different assay to measure 
the protein concentration, e.g. Bradford.  Alternatively, I could just trust my 
gel and see what happens when I add 10x more protein to the reaction.
