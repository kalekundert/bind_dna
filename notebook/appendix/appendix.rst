********
Appendix
********

Azure Sapphire optics
=====================
The Sapphire has slightly different lasers than the typhoon, but information 
about them was easily available in the control software.  Each laser has only 
one possible filter (this may actually be the case with the Typhoon, too, I'm 
not sure):

.. datatable:: sapphire_lasers_filters.xlsx

   I compiled the list of fluorophores for each laser by looking at page 12 of 
   the Sapphire brochure: :download:`20190723_confirm_cis_display_with_fluorescent_protein/sapphire_brochure.pdf`.  
   I also compared the imager's lasers and filters with the fluorophore's 
   excitation and emission spectra.  SYPRO Ruby doesn't seem like it'd be 
   excited by a 658 nm laser at all, but if Azure says it'll work, I'll trust 
   them.

Gel Doc EZ optics
=================
UV tray: 280-400 nm, max 300 nm
Blue tray: 430-460 nm, max 440 nm
White tray: visible spectra

Zif268 binding buffers
======================
Many different Zif268 binding buffers have been described.  Below are some of 
the ones I've found:

.. datatable:: zif268_buffers.xlsx

Some general takeaways:

- The buffer is typically ≈10 mM and around pH=7.5.
- The monovalent salt concentration is often ≈100 mM.  The specific cation 
  doesn't seem important.
- Mg is often (but not always) present at 1 mM.
- The Zn counter-ion doesn't seem important.
- The reducing agent is probably important.  The only Cys residues in Zif268 
  are involved in coordinating the Zn.  Obviously the protein won't work if 
  these residues form disulfide bonds, but too much reducing agent might also 
  interfere with Zn binding (I'm not sure exactly what oxidation state the 
  Zn/Cys need to have).
- A variety of nonspecific binding inhibitors are used.

I used variants of the [Lam2011]_ buffer for my initial experiments, because it 
was the only reference I found to express Zif268 with PURExpress.

PURE reaction components
========================
I couldn't find the composition of NEB's PURExpress buffer, but the buffer used 
in the original PURE system [Shimizu2001] is (50 μL for mass/quantity units):

- 9 mM magnesium acetate
- 5 mM potassium phosphate, pH 7.3
- 95 mM potassium glutamate
- 5 mM ammonium chloride
- 0.5 mM calcium chloride
- 1 mM spermidine
- 8 mM putrescine
- 1 mM dithiothreitol (DTT)
- 2 mM each ATP and GTP
- 1 mM each of CTP and UTP
- 10 mM creatine phosphate
- 2.8 A260 units tRNA mix (KBK: 112 ng/µL, 5.6 µg) (Roche, Mannheim, Germany)
- 0.5 μg 10-formyl-5,6,7,8-tetrahydrofolic acid
- 0.1 mM each of amino acids
- 12 pmol (32.4 μg, 240 nM) ribosome
- 1 μg IF1
- 2 μg IF2
- 0.75 μg IF3
- 1 μg EF-G
- 2 μg EF-Tu
- 1 μg EF-Ts
- 0.5 μg RF1
- 0.5 μg RF3
- 0.5 μg RRF
- 30–300 units of each ARS and MTF
- 0.2 μg creatine kinase (Roche)
- 0.15 μg myokinase (Sigma, St. Louis, MO)
- 0.054 μg nucleoside-diphosphate kinase
- 0.1 units pyrophosphatase (Sigma)
- 0.5 μg T7 RNA polymerase (NEB: 99 kDa; KBK: 100 nM)

PURExpress FAQs
===============
7. "The concentration of ribosomes in a standard reaction is approximately 2 µM 
   ±20%." (I calculated 2.4 µM in :expt:`30`.)

11. "For the control template DHFR, we estimate the ribosome recycled 5 times 
    successfully."

14. A T7 terminator is recommended even for linear templates, because it helps 
    improve mRNA stability.

Isoelectric points
==================
.. datatable:: zif_pi.xlsx

   Predicted isoelectric points (pI) from 
   `https://web.expasy.org/compute_pi/`__.

Nucleic acid stains
===================
I decided to use GelGreen for all applications:

GelGreen vs GelRed:

- Both dyes are excited by UV and blue light, but GelRed is much more excited 
  by UV and GelGreen is much more excited by blue light.  Since I will be using 
  the laser scanner when I really care, and the laser scanner only has blue 
  light, GelGreen makes more sense.

- GelRed is more sensitive for ssDNA and RNA than GelGreen is, but I think the 
  aforementioned advantages of GelGreen will make it the better stain in these 
  cases too.

GelGreen vs. SYBR

- I can't find any comparisons, but I think SYBR Gold is probably the most 
  sensitive dye in general, and the SYBR Green II is the most sensitive RNA 
  dye.  All of the SYBR dyes are well-excited by blue light.

- But the SYBR dyes seem much more difficult to work with: they need to be 
  stored at −20°C (and even then expire after a year), and they are sensitive 
  to the pH of the staining buffer (it needs to be pH 8).  They also aren't 
  considered "safe".

- GelGreen is also much cheaper.

FluoroTect GreenLys
===================
I want to know the length and molecular weight of the FluoroTect GreenLys 
reagent, for the purpose of properly labeling my gels.  Unfortunately, this 
information doesn't seem to be readily available, so I had to piece it 
together.

The `FluoroTect GreenLys manual`__ states that the reagent is based on the E.  
coli lysine tRNA.  There are `at least 3 such tRNAs`__: lysY, lysZ, and lyzQ.  
They all seem quite similar, so I chose to use lysQ purely because it showed up 
first in the search engine.  lysQ has:

- Length: 76 nt
- MW (tRNA only): 24.361 kDa 
- MW (with lysine and BODIPY): 24.769
  - BODIPY: 262 Da
  - lysine: 146

__ https://www.promega.com/-/media/files/resources/protocols/technical-bulletins/0/fluorotect-greenlys-in-vitro-translation-labeling-system-protocol.pdf
__ https://biocyc.org/ECOLI/NEW-IMAGE?type=LOCUS-POSITION&object=G6392&chromosome=COLI-K12&orgids=ECOLI&bp-range=780646/782584
