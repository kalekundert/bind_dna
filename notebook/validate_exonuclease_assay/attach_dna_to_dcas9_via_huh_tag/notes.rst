*******************************
Attach DNA to dCas9 via HUH-tag
*******************************

.. toctree::

  first_try/notes
  try_manganese/notes
  try_to_make_dna_homogeneous/notes
  confirm_dcas9_pcv2_concentration/notes
  try_assuming_2µM_dcas9_pcv2/notes
  try_longer_ori/notes
  try_fresh_aliquot/notes

[Lovendahl2017]_ describes a system for labeling proteins with nucleic acids.  
This system uses natural HUH-domains that react to form covalent bonds with the 
backbone of specific ssDNA sequences.  Although this method couldn't be used to 
make libraries, it can be used to validate the underlying exonuclease assay.

The first step is to show that I can efficiently attach DNA to a protein.  
Ultimately I want that protein to be Zif268, the model system I've been using, 
but I'm going to start with Cas9 because Jorge Marchand and others in the lab 
have already purified dCas9-PCV2 (PCV2 being a HUH domain).

Considerations
==============

Buffer
------
The buffer from [Lovendahl2017]_ is:

- 50 mM HEPES, pH=8.0
- 50 mM NaCl
- 1 mM MgCl₂
- 1 mM MnCl₂

This is slightly different from the buffer used for PCV2 in [VegaRocha2007]_:

- 20 mM Tris pH=7.4
- 100 mM NaCl
- 2.5 mM MnCl₂

Jorge says that he just uses CutSmart buffer, which is:

- 20 mM Tris, pH 7.9 at 25°C
- 50 mM KOAc
- 10 mM Mg(OAc)₂
- 100 µg/mL BSA

I also want a buffer compatible with the DNA-binding protein I'll be testing 
(e.g. dCas9 or Zif268).  NEB recommends the following buffer for Cas9:

- 20 mM HEPES, pH 6.5 at 25°C
- 100 mM NaCl
- 5 mM MgCl₂
- 0.1 mM EDTA

I would expect the [VegaRocha2007]_ buffer to work the best:

- Reasonable pH and salt concentration.
- Contains the best divalent metal ion: Mn²⁺
- No BSA: [Lovendahl2017]_ show that BSA is not needed.
- No EDTA: This would reduce the concentration of divalent metal.

.. update:: 2020/11/15

  Tris is known to chelate divalent metal cations, so in fact it might be 
  better to uses HEPES.

PBS supplemented with MnCl₂ might also work well.

I'm just going to use CutSmart for my initial experiments, because it 
apparently works for Jorge.


