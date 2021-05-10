************
Try PUREfrex
************

I learned in :expt:`67` that using PUREfrex instead of PURExpress stops the 
linker from being cleaved from the mRNA.  With that problem solved, I now want 
to test if I can observe attachment of mWasabi to its mRNA when using PUREfrex.

Try PUREfrex2.0 --- 2021/03/01
==============================
.. protocol:: 20210301_attach_via_purefrex.txt

.. figure:: 20210301_attach_mwasabi_purefrex.svg

Observations:

- I added 2.2x less RNase inhibitor than I should have, due to a bug in my 
  protocol script.  I don't think this should have much of an effect, though.

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

Try PUREfrex1.0 --- 2021/04/16
==============================
.. protocol:: 20210416_attach_via_purefrex.txt

.. figure:: 20210416_attach_via_purefrex1.svg

.. protocol:: 20210420_gel_laser_scanner.txt

.. figure:: 20210420_attach_via_purefrex1_urea.svg

  The control lanes are faint because I had <1 µL left over.

Observations:

- The SDS PAGE gel didn't run well, and the TBE/urea gel is also somewhat 
  distorted lower down.  I suspect that this has to do with with fact that (i) 
  I had samples with high salt concentrations and/or (ii) I had samples with 
  very different salt concentrations.

- The precipitation of the SDS PAGE loading buffer causes sample to be lost.  
  In the SDS PAGE gel, there is a clear loss of full-length mRNA at higher salt 
  concentrations.  In the TBE/urea gel, in contrast, every sample has the same 
  high level of full-length mRNA.  This indicates that full length mRNA is 
  getting trapped in the SDS precipitate.  Note that I saw the same effect in 
  :expt:`99` (2021/03/02).

  If I want to use SDS PAGE going forward, I'll need to desalt my samples.  
  Maybe I can get away without this for the FLAG peptide, since it's much 
  smaller and [Reyes2021]_ doesn't mention any desalting steps.  But in general 
  this is something I'll have to work around.

  .. update:: 2021/05/06

    I looked more closely into ordering desalting columns, and I realized that 
    they have molecular weight cutoffs.  Specifically, the Zeba micro spin 
    desalting columns claim to retain any molecules <1 kDa and elute any 
    molecules >7 kDa.  This means that I probably can't desalt my FLAG (2 kDa) 
    samples.  Zif268 (11 kDa) is on the border, but should be eluted.  mWasabi 
    should be totally fine.

- The TBE/urea gel might show the mRNA/protein fusion.  

  The band above the full-length mRNA is about what I'd expect the fusion to 
  look like, especially since the band grows more pronounced at high salt 
  concentrations.  There is no corresponding green band, but the green channel 
  is pretty much missing completely in this gel.  I think this may be some 
  consequence of the over-the-weekend −20°C incubation.

  I can't really quantify the intensities of the potential fusion bands, 
  because they're just too close to the full-length band.  But it does seem 
  qualitatively as if 375 mM KCl, 32.5 mM MgOAc condition gives the most 
  coupling.  This is the same condition that was optimal for [Reyes2021]_.  
  Note that I can't even estimate the fraction of coupled mRNA from this image, 
  because the full-length band is highly oversaturated.
  
- This gel is very similar to the PUREfrex2.0 gel.  But comparing the two gels 
  is still quite informative:

  - The level of mWasabi expression is much lower in PUREfrex1.0.  This is 
    consistent with what I saw in :expt:`99`, and consistent with the marketing 
    for PUREfrex2.0.

  - The same two GFP bands are present in both reactions.  This is noteworthy 
    because PUREfrex1.0—unlike PUREfrex2.0—only exhibits a single band when 
    expressing mWasabi from regular mRNA (:expt:`99`). With PUREfrex2.0, I just 
    attributed the second band to whatever mysterious process was creating an 
    extra band with regular mRNA.  With PUREfrex1.0, that explanation doesn't 
    work.  Furthermore, based on the fact that the MWs of the PUREfrex2.0 bands 
    didn't really agree, I suspect that something complicated is happening.

    - Is the 40 kDa band the mRNA display product?

      - What else could it be?

      - It's not the right MW.  I'd expect the mWasabi/mRNA fusion to run 
        slower than both the mRNA (≈270 kDa) and mWasabi (27 kDa), but this 
        band runs much faster than the mRNA.

        One possible explanation for this could be that DNA both has a smaller 
        hydrodynamic radius and a smaller charge than SDS-coated protein.  This 
        would allow DNA to migrate faster when attached to protein.  This idea 
        seems at least plausible given some estimates of hydrodynamic radius 
        (assuming both molecules adopt extended linear conformations):

        - SDS-coated protein: 18Å [`source 
          <http://hackert.cm.utexas.edu/courses/ch370/old2008/Electrophor/Electrophoresis.htm>`__.  
          This isn't a primary source, and gives the above number without any 
          reference.  So takes this with an extra-big grain of salt.

        - dsDNA: 11-13Å [`source 
          <https://bionumbers.hms.harvard.edu/bionumber.aspx?id=105243&ver=4&trm=dna+radius&org=>`__].  
          This doesn't account for possible secondary structure in the mRNA, 
          which would not be denatured in SDS-PAGE conditions.

        The relative charge of the protein vs. mRNA boils down to the density 
        with which SDS coats the protein vs. the intrinsic charge density of 
        the ssRNA backbone.  Assuming that SDS forms a micelle around the 
        protein, it's at least plausible that the proteins would have a higher 
        charge density.
      
      - The Cy5 signal doesn't really superimpose with the mWasabi signal.  
        This is more clear in the PUREfrex2.0 gel, which is just a higher 
        quality gel.  Even in the PUREfrex1.0 gel, the brightness of the 
        almost-superimposed red band is inversely proportional to the 
        brightness of the green band.

Conclusions:

- SDS precipitation is a major problem.  I'll have to either continue using 
  TBE/urea PAGE, or start desalting my samples.

- I might have seen coupling.  I'll have to repeat the experiment to see if I 
  can visualize the protein as well as the mRNA.

Without SDS PAGE --- 2021/05/03
===============================
.. protocol:: 20210503_attach_via_purefrex.txt

.. figure:: 20210503_attach_via_purefrex.svg

Observations:

- I do not see any coupling between the mRNA and the protein.
  
  - The PUREfrex 1.0 reactions have the same shifted bands that I saw in the 
    previous reaction (albeit more faint), but they are not superimposed with 
    any GFP signal.

  - The PUREfrex 2.0 reactions don't really have the same shifted bands, 
    although they do seem to have more smeared signal above the main mRNA band 
    in this high-salt conditions.

- All of the PUREfrex 2.0 reaction have a significant amount of Cy5 stuck in 
  the well.  I assume this is mRNA that's still associated with the ribosome, 
  but I don't know.

- I don't know which green bands are GP and which are FluoroTect.  In 
  retrospect, I probably shouldn't have used FluoroTect without having a 
  +FluoroTect −mRNA control.

- I don't know why the green bands don't exactly correspond between the 
  PUREfrex 1.0 and 2.0 conditions.

- The concentration of free linker in this f89 prep is surprisingly high.  
  Maybe this is because I decreased the concentration of ligase.  I should've 
  been more careful about that.

- The crystal violet loading dye worked really well.

Conclusions:

- These results are not promising.
