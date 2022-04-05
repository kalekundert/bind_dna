**********************
Measure Zif268 binding
**********************

2022/03/17
==========

.. protocol:: 20220310_b1h.pdf 20220310_b1h.txt

.. figure:: 20220317_b1h.svg

Observations:

- The control didn't work well this time, compared to :expt:`41`.  
  Specifically, the negative control (AAA target) grew quite well in the 
  selective condition.  The cells also grew 10-20x better in the permissive 
  condition than they had previously.

  - I wonder if there was cross-contamination between my cultures, e.g.  
    because I was using a 24-well block instead of culture tubes.

  - I could sequence the colonies to find out: probably worth doing.

- The single plasmid construct had a range of activities.  Most active was 
  s16/s22, which had the insert before AmpR and no barcode.

- I only had +/− barcode constructs for the "before AmpR insert" prepared for 
  this experiment.  Surprisingly, the barcode made a difference: Without it the 
  selection was extremely effective; with it less so.

  - One explanation is that this construct has higher AmpR expression.  The 
    barcode insert includes bidirectional terminators, which in this case 
    prevent the AmpR gene from being included in the Zif268-rpoZ transcript.  
    The strain without this insert should have higher AmpR expression, which 
    could explain why it seems healthier.  See :expt:`161` for more on this 
    idea.

- I can't really understand why inserting the Zif268-rpoZ gene into the f1 ori 
  would be substantially better than inserting it between the f1 ori and the 
  pSC101 ORI.  I want to sequence and/or repeat things before I read too much 
  into this.

Next steps:

- Sequence colonies
- Repeat with all constructs: are results reproducible?

2022/03/23
==========
Today I'm going to repeat the experiment with only the controls and s16/s22.  
My hope is that doing the experiment with fewer samples (and using individual 
culture tubes) will give better results.

.. protocol:: 20220323_plate_assay.pdf 20220323_plate_assay.txt 20220325_b1h_image_processing.txt

.. figure:: 20220325_b1h_s4_s5_s16_s22.svg

Observations:

- The controls worked well this time.  I suspect the problem with the previous 
  experiment was contamination.

- s16 actually performed better than the controls: you can see that the 
  colonies are bigger (although still not quite as big as the positive 
  control).  This could be a sign that even in the original constructs, the 
  antibiotic resistance gene is not being expressed highly enough.

- Reverse pipetting (which I used on the s16/s22 selective plate) completely 
  eliminated the "splash colonies" that are visible on every other plate.

- I did the image analysis differently, most significantly (i) subtracting 
  background signal and (ii) setting contrast based on the number of saturated 
  pixels.  I think this makes it a lot easier to see/compare the colonies.

Nest steps:

- I'm hesitant to repeat this assay with lots of samples again, because it was 
  difficult and I might have contamination again.  But I want to retest the 
  constructs with Zif268 inserted after/within f1.  I'll either do that in 
  smaller batches with this assay, or perhaps all at once with the plate reader 
  assay (once I've analyzed the data and decided whether or not to trust it).

2022/04/01
==========

.. protocol:: 20220329_plate_assay.txt 20220325_b1h_image_processing.txt

.. figure:: 20220401_b1h.svg

Observations:

- The s14/s20 results are consistent with what I saw in the 3/17 experiment.  
  The controls were wrong in that experiment, perhaps because of 
  cross-contamination, but maybe the rest of the data are correct.

- Including the barcode and the terminators (s11/s17) seems to really break the 
  assay.  The strain actually grows very well in the selective condition with 
  the off-target site, and actually grows worse with the target site.  I saw 
  the same thing in the 3/17 experiment for the s12/s18 plasmids, and although 
  I don't necessarily trust those results, they match these results almost 
  perfectly.  These two pairs of plasmids are very similar, with the former 
  having the insert within the f1 region and the latter having it after.  So it 
  seems likely that something about having the barcode/terminator shortly after 
  the HIS/URA genes is really bad.  Some ideas:

  - Cryptic promoters?  There are some predicted promoters in the zif-rpoZ 
    insert, including one that overlaps the L3S2P21/barcode (and so would only 
    be present in the barcode/terminator constructs).  But it's not 
    particularly strong, and the polymerase would have to go all the way around 
    the plasmid (including through several terminators) to transcribe the 
    HIS/URA gene.  I don't think this is relevant.

  - Increased mRNA stability?  Do strong terminators give rise to more stable 
    mRNA, since there would be a stronger hairpin on the end?  This might be 
    true [Ahn2008]_.  I think this might be why these strains grow so well in 
    the selective condition.
    
    This explanation doesn't work as well for the "after f1" constructs 
    (s12/s18), though, since they still have the rrnB T1 terminator after the 
    HIS/URA genes.

  - Antibiotic resistance?  If some AmpR expression is being driven by a 
    cryptic promoter upstream of the target site, it could be that Zif268 
    binding its target site would actually reduce AmpR expression.  I think 
    this might be the reason why the target strain grows worse than the 
    non-target strain.
  
Next steps:

- Finish testing all these strains.

- Think about the terminator/mRNA stability of the HIS/URA transcript.

  - What's happening in the original assay?

  - Include a ribozyme just before the terminator?

  - qPCR for HIS/URA gene?

    - Might be the best way to answer this question, and it tells me if the 
      RNAseq assay is a viable idea.  That aid, it should be the case already 
      that survival is correlated with mRNA expression.  This wouldn't directly 
      tell me *why* there is more mRNA.

    - I've done this assay already (for the sgRNA project).  This basic 
      protocol is: TRIzol, RT (with random hexamers), qPCR.

2022/04/02
==========

.. protocol:: 20220330_plate_assay.txt

.. figure:: 20220402_b1h.svg

Observations:

- This result is not consistent with the 3/17 experiment, in which the s18 was 
  viable even in selective conditions.  I suspect that the 3/17 experiment had 
  contamination.  The current result is also more consistent with my terminator 
  hypothesis.

Next steps:

- Test s15/s21

- Maybe test s12/s18 and s11/s17 together.

- Try best AmpR promoter with p194.


