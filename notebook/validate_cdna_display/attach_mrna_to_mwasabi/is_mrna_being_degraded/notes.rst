***********************
Is mRNA being degraded?
***********************

2020/09/09:

In all of my reactions so far, I've seen that the linker appears to go from 
being fully attached to the mRNA (before the PURExpress reaction) to fully 
separated (after the ligation).  I can think of two broad explanations for 
this:

1. The mRNA is being completely degraded, leaving only the DNA linker.  The 
   problem with this explanation is that the degradation would have to be 
   extremely efficient, as there is no smearing evident.  It's hard to explain 
   how this could happen: The reaction has RNase-inhibitor, and shouldn't have 
   any RNase in the first place.

2. The linker is not actually ligated to the mRNA.  Either it was never ligated 
   in the first place, or something is undoing the ligation.  The problem with 
   this explanation is that it's just implausible.  Experiments such as 
   :expt:`50` convincingly show that the ligation is working,and I just can't 
   imagine how the ligation could be reversed.

.. update:: 2020/10/27

  George proposed a reasonable mechanism for the second possibility: RNaseH 
  contamination.  This would cleave the linker from the mRNA.

I think it's an important first step to learn whether or not the mRNA is 
actually being degraded. 

Considerations
==============

Methods to visualize RNA
------------------------
My initial thought was to try visualizing the mRNA with direct GelGreen 
staining.  However, after looking back at :expt:`18`, I don't think this is 
going to work.  There's a big smear in PURExpress that completely obscures the 
mRNA band.  Some ideas:

- Use a urea gel.

  - The gel in :expt:`18` was an SDS PAGE gel, and it's probably true that a 
    urea gel would resolve more bands.  But this seems like a slim hope.

- Add the spinach/broccoli aptamer to the mRNA [Filonov2015]_.

  - According to [Filonov2015]_, the spinach and broccoli aptamers can be 
    detected in both native and denaturing gels.  A brief wash step is used in 
    either case.

  - I'll have to put the aptamer on the 5' end of the mRNA, before the RBS.  
    [Filonov2015]_ claims that the aptamer is not particularly sensitive to its 
    context, so there's nothing a priori wrong with putting it on the 5' end.  
    The 5' end is easier for cloning too.

  - If the mRNA is being degraded, the aptamer will also be degraded and I 
    won't be able to see it.  But I should still be able to distinguish 
    degraded from not-ligated mRNA.

  - The reason to prefer broccoli over spinach is that the former is shorter 
    (47 nt).  Both seem to perform equally well.

  Nomenclature:

  - DFHBI-1T: Di Fluoro Hydroxy Benzylidene Imidazole 1-Trifluoro (That leaves 
    out a whole bunch of words, but explains what the letters stand for)

  Protocol:

  - Run gel

  - Wash 3x 5 min in water.

  - Stain for 10-30 min in following buffer:

    - 10 µM DFHBI-1T
    - 40 mM HEPES (pH 7.4)  (pH>6 is `important 
      <https://en.wikipedia.org/wiki/Spinach_aptamer>`__)
    - 100 mM KCl
    - 1 mM MgCl₂

  Prices:

  - DFHBI: $185 per 10 mg
  - DFHBI-1T (2x brighter): $268 per 10 mg

- Northern blotting

  - I already ordered most of the things I would need, including the probe 
    (o135 targets mWasabi).

  - One problem is that the membranes used for blotting autofluoresce at the 
    same wavelengths as Cy5.  This means that I'd probably end up having to 
    visualize the gel before blotting to see where the linker (and mWasabi) 
    are.

After talking with Fitzy about it, I think I'm going to try to pursue the 
Northern blotting and broccoli aptamer ideas in parallel.

Broccoli vs. tBroccoli
----------------------
[Filonov2015]_ uses both broccoli and tBroccoli (broccoli fused to a tRNA 
scaffold to help improve folding).  I decided to use plain broccoli because 
[Filonov2015]_ shows that the tRNA scaffold is recognized by certain RNases in 
bacteria (specifically RNases E, T, and PH).  The broccoli aptamer seems to 
fold just as robustly, and it's shorter.

Sequence: AGACGGTCGGGTCCAGATATTCGTATCTGTCGAGTAGAGTGTGGGCT


Results
=======

Visualize mRNA via broccoli --- 2020/11/09
------------------------------------------

.. figure:: 20201109_confirm_dfhbi_stain.svg

.. datatable:: 20201109_confirm_dfhbi_stain.xlsx

- The DFHBI-1T staining worked.

- I got much lower sensitivity than [Filonov2015]_.

  The smallest quantity of mRNA I could detect on this gel was 8 ng.  
  [Filonov2015]_ claims that quantities as low as 100 pg can be detected.  This 
  implies that my signal is about 80x less than it should be.

  One likely problem is that I used Tris instead of HEPES for the staining 
  buffer.  (I just ran out of time to make the HEPES buffer.)  Tris is known to 
  chelate magnesium and many other metal ions [Fischer1979]_.  A quote from the 
  abstract: "great reservations should be exercised in employing Tris as a 
  buffer in systems which also contain metal ions".  This may have had the 
  effect of reducing the amount of magnesium available to the broccoli aptamer, 
  preventing it from folding correctly.
