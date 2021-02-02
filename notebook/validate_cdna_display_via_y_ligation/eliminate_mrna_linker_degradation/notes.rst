*********************************
Eliminate mRNA-linker degradation
*********************************

.. toctree::
   :hidden:

   check_degradation_via_broccoli/notes
   detect_rnase_h_activity/notes

2020/09/09:

In all of my protein expression reactions with ligated mRNA-linker so far 
(:expt:`65`, :expt:`61`), I've seen that the linker appears to go from being 
fully attached to the mRNA (before the reaction) to fully separated (after the 
reaction).  I can think of two broad explanations for this:

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
  contamination.  This would cleave the linker from the mRNA.  See :expt:`77`

The goal of this series of experiments is to diagnose and solve this problem.

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



