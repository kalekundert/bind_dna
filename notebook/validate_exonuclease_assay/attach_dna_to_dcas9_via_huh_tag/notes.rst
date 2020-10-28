*******************************
Attach DNA to dCas9 via HUH-tag
*******************************

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

PBS supplemented with MnCl₂ might also work well.

I'm just going to use CutSmart for my initial experiments, because it 
apparently works for Jorge.


Methods
=======

Prepare DNA targets
-------------------
.. protocol:: 20200225_make_dna_beads.txt 

.. figure:: 20200225_huh_tag_targets.svg 

- The reactions were very clean.

- Yields were lower than I expected, 20 µL at ≈50 ng/µL.  The product is 
  relatively small; I wonder if PCR cleanup would give better yield than the 
  beads (or if I should've used more beads).

Attach DNA --- 2020/09/18
-------------------------
.. update:: 2020/10/15

  I analyzed this data thinking that the bright protein band was dCas9, when it 
  is realy BSA.  The fainter band on top is dCas9.

.. protocol:: 20200918_attach_huh.txt

.. figure:: 20200918_attach_huh_dna_to_dcas9.svg

- The DNA is very faint.  I think the problem is that I didn't rinse the gel 
  before staining it.  I noticed that the staining solution I used now foams up 
  when I shake it, suggesting that it picked up some SDS from the gel.  I also 
  noticed that the staining solution stayed more orange than usual after 
  subsequent uses, suggesting that not as much GelGreen was absorbed by the gel 
  (presumably because the SDS was preventing it from binding DNA).

- Assuming that the DNA is the band that changes between each condition (and 
  that the sgRNA is the band that stays in place), it's bizarre that the DNA 
  shifts, but doesn't run exactly with dCas9 in either condition.

- I'd like to see how the DNA runs by itself.  As it is, it's kinda hard to 
  distinguish the DNA from the sgRNA.  I could also try the reaction without 
  the sgRNA; the HUH-domain should still react.

- I'd also like to see a DNA ladder (and a protein ladder, for that matter).

Attach DNA --- 2020/10/13
-------------------------
.. update:: 2020/10/15

  I analyzed this data thinking that the bright protein band was dCas9, when it 
  is realy BSA.  The fainter band on top is dCas9.

.. protocol:: 20201013_make_custom.txt 20201013_attach_huh.txt

.. figure:: 20201013_attach_huh_dcas9.svg

- The DNA channel has the expected shifts, but the dCas9 channel doesn't have 
  any shifts at all.  I really don't get how that's possible.  How can the DNA 
  shift without the protein?  Some thoughts:

  - I should check to make sure that the previous refs did SDS PAGE...  Update: 
    They did.  And they only stained with Coomassie (i.e. think didn't stain 
    the DNA), so I definitely should see something in the Coomassie channel.

  - dCas9 is 160 kDa and f16 is 255 kDa, so it's definitely reasonable to think 
    that coupling these two molecules would have some effect.

  - Maybe what I'm seeing is binding without covalent attachment.  I saw a 
    similar effect with my :expt:`Zif268 EMSA experiments <35>` experiments, 
    i.e. the DNA was shifted behind where the protein normally runs, but the 
    protein itself was not shifted.  That was a native gel, though, and this is 
    an SDS gel: any non-covalent binding should have been eliminated.

    If this is somehow the case, it may help to actually use manganese in my 
    reaction buffer.

- I noticed that I may have ordered the wrong primer sequence:

  ======================================================  ========================================================
  Source                                                  Sequence
  ======================================================  ========================================================
  o102                                                    ``TAAAGTATTACCAG/iSp9/atctttctgacgcagatgaa``
  Jorge (via email)                                       ``TAAAGTATTACCAG(NNNNNNNNNNNNNNN)/iSp9/PRIMER``
  [Lovendahl2017]_, Table 1                               ``aagtattaccagaaa``
  [Lovendahl2017]_, Table S18, Donor-quencher oligos      ``IowaBlackFQ/AAAGTATTACCAGA/FAM``
  [Lovendahl2017]_, Table S18, Amino oligos               ``AAGTATTACCAGAAA/NH2``
  [VegaRocha2007]_, P10                                   ``AAGTATTACC``
  [VegaRocha2007]_, P12                                   ``AAGTATTACCAG``
  ======================================================  ========================================================
 
  Looking at all these sequences, I can see why I'd expect o102 to work: both 
  [VegaRocha2007]_ sequences are even shorter!  But I can also see why this 
  might cause problems: all of the [Lovendahl2017]_ sequences are longer, and 
  Jorge recommended a bunch of spacer nucleotides!  It's also not ridiculous to 
  think that the fusion might need a little more space than the free protein 
  for some reason.

  I was thinking about ordering primers for a different length amplicon anwyays 
  (either shorter or longer, see below), so maybe I should just do that with 
  Jorge's primer.

- Note that f12 runs a bit slower than f16.  This is presumably because f12 is 
  15 nt longer.

- dCas9 appears to be faintly visible in the GelGreen channel, even without 
  sgRNA.  I suppose it's not unreasonable to think that GelGreen could 
  intercalate somewhere in a big protein like Cas9.

- It's unfortunate that free DNA seems to run pretty much the same as free 
  Cas9.  I wonder if I should try this experiment with a much shorter or longer 
  amplicon.  Because of how the primers are designed, I'd have to change the 
  forward primer to change the size of the amplicon.  The resulting DNA would 
  not be ideal for the actual protection assay, but would be better for simply 
  confirming attachment.

  - With the primers I have on hand, I can only get about as short as 300 bp.

  - I can get 1.5 kb by using o112.  Otherwise my options are 550 bp or 2.2 kb.

  - Maybe 300 bp and run for another 20 min or so...

- I could run pure Cas9 as a control.  I'd have to buy it though.

