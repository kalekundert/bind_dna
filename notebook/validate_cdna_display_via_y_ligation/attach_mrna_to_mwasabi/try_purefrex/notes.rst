************
Try PUREfrex
************

I learned in :expt:`67` that using PUREfrex instead of PURExpress stops the 
linker from being cleaved from the mRNA.  With that problem solved, I now want 
to test if I can observe attachment of mWasabi to its mRNA when using PUREfrex.

Results
========
.. protocol:: 20210301_attach_via_purefrex.txt

.. figure:: 20210301_attach_mwasabi_purefrex.svg

Observations:

- The image is really poor quality.  I don't know why.  I tried restarting the 
  software and the imager, but neither had any effect.

- None of the mWasabi appears to have reacted with the mRNA.

  There are two bands in the mWasabi channel.  I saw the same two bands in 
  :expt:`67`.  The gels in :expt:`99` (mWasabi with stop codon) also have two 
  bands.  The lower bands in both experiments seem to be about the same size 
  (≈30 kDa) but the upper bands are notably different (≈40 kDa vs ≈32 kDa).

  The upper band is more faint at high salt concentrations.  This could be due 
  to the precipitate, as only the highest salt concentration had a significant 
  amount of precipitation.

- There is still a significant amount of RNA degradation in the PUREfrex 
  reactions.  The −PUREfrex control, in contrast, looks pretty clean.

  The degradation seems to be accelerated at higher salt concentrations.  I can 
  think of two explantions for this:
  
  - RNase A is somehow activated by K⁺ or Mg²⁺.
  - The salts themselves are contaminated with nuclease.
    
  Some follow-up experiments I could do:
  
  - Incubate f89 alone with salt.  If the nuclease is in the PUREfrex, I should 
    observe no degradation.  I could also setup separate KCl and MgOAc 
    incubations, to see if the effect is due to one salt in particular.

  Notably, [Reyes2021]_ state that they observe no RNA degradation.  It's 
  possible that PUREfrex1.0 has less nuclease contamination than PUREfrex2.0.

- All of the +PUREfrex conditions have a clear Cy5 band just above the free 
  linker band.  I wonder if this the maximally degraded mRNA product (which 
  runs slower than free linker because it still has a little bit of 
  RNA---either because the nucleases only cut at certain sequences, or because 
  the DNA/RNA duplex is protected from nuclease activity.

- [Reyes2021]_ used MgCl₂, not MgOAc.  It probably doesn't make a big 
  difference, but it's something to be aware of.  I checked the various salt 
  concentrations I compiled in :expt:`62`, and both MgCl₂ and MgOAc have been 
  used.  [Naimudden2016]_ used MgCl₂ as well, though.

- The gel looks pretty much identical after the −20°C incubation.

Check salt degradation --- 2021/03/02
=====================================
Based on the results from yesterday, I wanted to see if the salt itself was 
responsible for the degradation.

.. protocol:: 20210302_check_salt_degradation.txt 20210302_gel.txt

.. figure:: 20210302_check_salt_degradation.svg

Observations:

- The MgOAc bands are missing in the SDS PAGE gel.  Note that I observed SDS 
  precipitate only in the KCl samples, as expected, so this is not attributable 
  to precipitation.  It's also not attributable to any kind of degradation, 
  because the Cy5 moiety is entirely gone.

  I don't know what happened, but since all the bands are present in the urea 
  gel, I'm going assume that this is some kind of artifact.

- The +salt lanes in the urea gel show no degradation relative to the −salt 
  lane, so I think it's fair to conclude that nothing in the salt is 
  responsible for the degradation I've seen in my PUREfrex reactions.  

  I don't really know what else I can do, though, other than trying PUREfrex1.0 
  and hoping it has less nuclease activity...


