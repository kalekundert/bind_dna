********************************
Release repA with rRNA digestion
********************************

2021/02/04:

The results from :expt:`36` are consistent with repA binding the oriR sequence 
as intended, but the protein doesn't appear to be freely soluble or monomeric.  
Instead, it seems to be associated with large PURExpress components.

Considering the particular hypothesis that repA is associated with the 
ribosome, I want to test if incubating the reaction with RNase frees the 
protein.  See :expt:`31` for an earlier discussion of trying to degrade the 
ribosome with RNase.  This experiment did not seem to be successful, but it was 
focused on purifying repA.  Electrophoresis is simpler than purification, so I 
think it's possible that I'll see something here that I didn't see previously.

Results
=======

EMSA --- 2021/02/04
-------------------
.. protocol:: 20210204_digest_ribosomes.txt

.. figure:: 20210204_release_with_rnase_37_65.svg

Observations:

- The lanes look really smeary and wavy.  One possible explanation is that the 
  concentration of glycerol in the samples is too high, on account of adding 
  the RNase cocktail.  It's a little hard to believe that 0.5 µL would have 
  that big of an effect, though.

  - I wonder if the sample have enough glycerol for me to just not use any 
    loading dye at all...

- I don't know what happened to the DNA in the −mWasabi,−oriR lanes.  I can't 
  think of any reason why the DNA would get stuck in the well in that 
  condition...

  I wonder if I somehow mixed up the f62 and f63 lanes?  I really don't think I 
  did, but the data would make more sense that way.

Conclusions:

- I can't really take away much from this gel.

- Maybe I should try the agarose again, on the theory that digesting the 
  ribosomes may allow things to move backwards.  I don't really think that 
  would give me good data, though.

- Maybe it wasn't great to exclude the RNase inhibitor.
