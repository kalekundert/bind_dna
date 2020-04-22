********
Appendix
********

Azure Sapphire optics
=====================
The Sapphire has slightly different lasers than the typhoon, but information 
about them was easily available in the control software.  Each laser has only 
one possible filter (this may actually be the case with the Typhoon, too, I'm 
not sure):

.. datatable:: 20190723_confirm_cis_display_with_fluorescent_protein/sapphire_lasers_filters.xlsx

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

PURE reaction components
========================
I couldn’t find the composition of NEB’s PURExpress buffer, but the buffer used 
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
- 12 pmol (32.4 μg) ribosome
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
- 0.5 μg T7 RNA polymerase

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
