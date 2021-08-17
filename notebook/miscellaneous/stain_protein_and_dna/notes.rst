*********************
Stain protein and DNA
*********************

I've had problems staining SDS-PAGE gels to visualize both protein and DNA 
bands.  For example, see :expt:`71`.  This would be very useful to be able to 
do reliably, so I want to spend a little bit of time trying to get it right.

There are a number of things I could try:

- Staining order
- With/without microwaving
- Protein stain: Coomassie, SafeStain, SYPRO Ruby, SYPRO Orange

Maybe the methanol/ethanol in the Coomassie staining solution is causing 
problems.

.. update:: 2020/11/15

  I realized today that the pH of the Coomassie staining solution is very 
  acidic: pH<1 for standard staining, and pH≈2.1 for SimplyBlue SafeStain.  
  I've found anecodotal sources claiming that the phosphate backbone is 
  hydrolyzed at pH<1.  With the solution also being microwaved to near boiling, 
  it's not hard to imagine that the DNA is being significantly hydrolyzed.

  Note that duplex DNA is also destabilised at pH<5, due to protonation of the 
  bases.  That could also affect GelGreen binding.

  If this is the case, the best course of action would be to use a different 
  protein stain:

  - SYPRO Ruby:

    - Still requires acetic acid fixation.

  - SYPRO Orange/Red:

    - Acetic acid improves sensitivity, but does not have to be used.
    - SYPRO Red uses the 532 nm laser.
      - Certainly crosstalk with the 488 nm (FITC) channel
      - Maybe crosstalk with the 658 nm (Cy5) channel
    - SYPRO Orange uses the 488 nm laser:
      - Unlikely to crosstalk with the 658 nm (Cy5) channel.
    - So neither stain is compatible with GelRed/GelGreen, but SYPRO Orange 
      might be compatible with my Cy5-labeled oligos.
    - For EMSA experiments (e.g. native PAGE), would need to soak the gel in 
      0.05% SDS prior to staining: 
      https://www.thermofisher.com/us/en/home/references/molecular-probes-the-handbook/protein-detection-and-proteomics-technology/detection-of-the-total-protein-profile-in-gels-on-blots-on-microarrays-and-in-capillary-electrophoresis.html

  - SYPRO Tangerine:

    - No fixation
    - Seems likely to have significant crosstalk with GelGreen.

  - No-stain:
    
    - No explicit fixation step, but the composition of the reaction buffer 
      isn't disclosed.
    - Would use green laser; crosstalk may be an issue.

GelGreen + Coomassie --- 2020/11/04
===================================
.. protocol:: 20201104_test_staining.txt

.. figure:: 20201104_test_staining_order.svg

- The wash steps help GelGreen bind DNA.

- I didn't add enough BSA.  I should go back to :expt:`46` and work out how 
  much was in the buffers I was using.

- Coomassie staining totally eliminates GelGreen staining, whether it happens 
  before or after.

  I'm slightly curious which component of the Coomassie staining is 
  responsible.  I don't actually know what's in SimplyBlue, so I can't really 
  answer this question.  But it's interesting thinking about the components of 
  a traditional Coomassie stain:

  - methanol: Maybe GelGreen is very soluble in methanol, so much so that it 
    won't bind DNA.  This would require that significant amounts of methanol 
    remain in the gel after Coomassie staining.  Maybe I could try doing a more 
    extensive microwave wash after Coomassie staining and before GelGreen 
    staining  (This time I just shook in water for 10 min, but maybe that's not 
    enough).

  - acetic acid: According to Biotium, GelGreen binds DNA via both 
    intercalation and electrostatic interaction.  By lowering the pH, acetic 
    acid could interfere with electrostatic binding.  Maybe for this it would 
    help to have the GelGreen in a buffered solution, e.g. PBS.

  - Coomassie: Maybe GelGreen is very efficiently quenched by the background 
    levels of Coomassie.

Cy5 + SYPRO Orange --- 2021/08/17
=================================
In :expt:`122`, I tried simultaneously imaging DNA and protein with Cy5 and 
SYPRO Orange.  The protein wasn't visible in that gel, but I think that's just 
because I didn't add enough protein (and SYPRO Orange never seems to stain my 
ladder for some reason).  More importantly, the DNA was clearly visible.  I 
wasn't sure if this would be the case, because I thought that the acetic acid 
in the staining solution might damage the DNA.  Given that success, I want to 
take a closer look at this visualization protocol to see how well it works.

It's worth noting that this protocol isn't exactly what I'm looking for, 
because Cy5 isn't a stain.  But labeling DNA with Cy5 is usually not too hard, 
so I can imagine this protocol being useful in plenty of circumstances.
