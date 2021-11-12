**************************
Try bacterial IVTT systems
**************************

2020/09/10:

So far, I've been using PURExpress for all of my cDNA display experiments.  
However, in light of the fact that PURExpress seems to be negatively affecting 
my mRNA (:expt:`66`), I think it's worth trying some other systems:

- PUREfrex: Another commercial PURE system sold by GeneFrontier.  According the 
  [Doerr2019]_ and [Niederholtmeyer2013]_, it has no detectable nuclease 
  activity (while PURExpress does have some).  It also has no His-tagged 
  components, so I can do His-tag purifications (but not reverse 
  purifications).

  Unfortunately, PUREfrex is not available in B2P, so ordering might be a 
  hassle.  I requested a quote.  See this page: 
  https://www.genefrontier.com/en/solutions/purefrex/lineup/

- Promega S30 extract, NEBExpress: Two commercial S30 lysates.  PURExpress 
  should have lower nuclease activity than any lysate, but maybe there's 
  something else going on.  Additionally, most mRNA display systems use lysates 
  (specifically, rabbit reticulocyte lysate), so it's not like lysates are 
  completely off the table.  Really I just have both on hand, so I feel like I 
  might as well try them.

I think it's worth trying rabbit reticulocyte lysate as well, but that would be 
a more involved experiment because I'd first have to reclone my constructs.

2020/02/22
==========
.. protocol:: 20210222_compare_ivtt.txt

.. figure:: 20210222_compare_ivtt.svg

Observations:

- My f89 stock is not concentrated enough to reach the ideal concentration for 
  the NEBExpress reaction (400 nM).  The most I can reach is:

  .. math::

    \frac{\pu{1000 nM} \times \pu{1.44 µL}}{\pu{6 µL}} = \pu{240 nM}

  According to the results from :expt:`98`, that should be enough to get ≈50% 
  of the maximum expression.  Hopefully that will be good enough.

- I mistakenly entered the wrong target concentration for the PUREfrex 
  reaction.  As a result, I used 200 nM and not 500 nM.  I couldn't have 
  actually reached 500 nM anyways, but I could've gotten to 350 nM.  200 nM 
  should still be enough to get 70% yield, though.

- Not all of the mRNA from the PURExpress reaction was degraded.  :expt:`61` 
  actually exhibits the same effect.  My guess is that if I were to do a 
  timecourse, I'd find that the mRNA has a half-life of ≈15 min in PURExpress.

- The PUREfrex reaction shows a significant amount of traditional RNA 
  degradation.  I wonder if adding an RNase inhibitor would help with this.

  Before doing this, I think it'd be smart to re-optimize the amount of mRNA to 
  use in this condition.  It stands to reason that adding an RNase inhibitor 
  would increase the effective concentration of the RNA, which could affect the 
  amount of RNA that should be added.  This is also consistent with the fact 
  that PUREfrex seems to require significantly more mRNA than PURExpress 
  (:expt:`18`, :expt:`99`).

  I saw in :expt:`99` that expressing mWasabi with PUREfrex gave two 
  fluorescent bands.  I wonder if adding RNase inhibitor would also eliminate 
  that behavior...

- I can't really identify any of the green bands, because I didn't include 
  −mRNA controls, and it's certainly possible (especially for the extracts) 
  that the IVTT reactions contain fluorescent components.  I should include 
  −mRNA controls next time I do something like this.

- The extracts (NEBExpress and Promega S30) seem to have significant levels of 
  ssDNA nuclease activity.
  
  The 17 kDa Cy5 band in the PURExpress reaction corresponds to the free 
  linker.  This is consistent with (i) my hypothesis that the linker is being 
  freed from the mRNA by RNase H activity, (ii) the linker-only (o129) control 
  from :expt:`95`, and (iii) the MW of the linker, which is 19.9 kDa.

  In the extract lanes, the Cy5 band has shifted to 3 kDa.  Given the rough 
  correspondence between the protein ladder and the molecular weights of the 
  DNA components, this seems like it must be smaller than even the individual 
  oligo containing Cy5 (o126, 34 nt, 11.2 kDa).  This would only be possible if 
  these reactions contain ssDNA nuclease activity.  If I wanted to look into 
  this more closely, I could use o126 (+/- nuclease, even) as a control.  That 
  said, these extracts do not seem like a promising avenue.

PUREfrex with inhibitor --- 2021/02/26
======================================
.. protocol:: 20210226_purefrex_with_inhibitor.txt

  Densiometry analysis:

  - Subtract background

  - Divide each lane into four sections:

    - Stuck in the well
    - Full length mRNA
    - mRNA degradation products
    - Free linker

    We're interested in the ratio between the 2nd and 3rd sections.  In ImageJ, 
    I drew a vertical line through all lanes, so that the divisions between 
    these sections are at the exact same pixel offsets for each lane.

.. figure:: 20210226_purefrex_with_inhibitor.svg

Observations:

- The amount of mRNA degradation is qualitatively (and quantitatively) less for 
  the +inhibitor reaction.  This is consistent with some level of RNase A/B/C 
  contamination.

- The −PUREfrex control is saturated.  I included it in the densiometry 
  calculations anyways, but I need to keep in mind that those numbers will be 
  underestimates.

- The +PUREfrex lanes have a significant amount if signal that doesn't enter 
  the gel.  I suspect that this is mRNA that's still bound to the ribosome 
  somehow.  It could also be mRNA that is actually attached to GFP, although 
  there's no signal in the green channel.  I didn't include this signal when 
  calculating the percentage of full-length mRNA, although it's reasonable to 
  think that full-length mRNA would be more likely to bind the ribosome (e.g.  
  mRNA without the RBS would be less likely to bind).  So the percentage of 
  full-length mRNA might be as high as 10% and 34% respectively for the −/+ 
  RNase inhibitor conditions.

- The +inhibitor condition has slightly higher protein expression, which makes 
  sense.  Note that there are 2 (faint) green fluorescent bands, neither of 
  which appear in the −mRNA reactions.  Both are included in the "% expression" 
  calculation.

Conclusions:

- When I eventually do the exonuclease assay, I wonder if I can find a way to 
  put the target sequence/barcode on the 3' end of the mRNA, where it should be 
  more likely to end up attached to the protein even if the mRNA is degraded to 
  some extent.

Discussion
==========
- PUREfrex seems to have no appreciable RNase H activity.

- When RNase inhibitor is added to the PUREfrex reaction, relatively little of 
  the mRNA is degraded.
