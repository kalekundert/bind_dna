***************************
Express mRNA via PURExpress
***************************

In :expt:`19` (2020/02/18), I saw poor expression of protein from mRNA.  Here I 
want to determine how much mRNA is needed to get appreciable expression with 
PURExpress.  I decided to start with mWasabi, because it should be pretty easy 
to detect its expression.  I even think it might remain fluorescent in after 
SDS-PAGE based on this `high school/undergrad lab kit 
<https://www.bio-rad.com/en-us/product/pglo-sds-page-extension?ID=a41608e9-b348-43e0-98bb-d0ae12664e06>`__, 
which would be extra convenient.  Even if that doesn't work, though, it's still 
a protein that I know to be well-expressed.  DHFR (the NEB control) is another 
option.
  
Considerations
==============
NEB recommends using 1--5 µg mRNA per 25 µL PURExpress reaction.  This 
corresponds to 1.6 µL of 10 µM f11 per 10 µL reaction.

Results
=======

Express mWasabi --- 2020/02/25
------------------------------
.. protocol:: 20200225_split_purex_page.txt

.. figure:: 20200225_express_mwasabi_mrna.svg

   (A) Direct imaging of mWasabi fluorescence. mRNA, DNA: Concentrations refer 
   to the "template mRNA" stock solution in step 2 of the above protocol.  
   PURExpress: Steps 2-3 in the above protocol.  (B) Coomassie staining, zoomed 
   in on mWasabi.  (C) GelGreen staining.

- The 1 µM reaction seems to have produced significantly more mWasabi than the 
  10 µM reaction.  This effect is clear both with mWasabi fluorescence and with 
  Coomassie staining, which means it isn't due solely to SDS interfering with 
  mWasabi fluorescence.  The may be some interference by SDS---relative to the 
  DNA template, the 1 µM mRNA template does seem brighter when imaged via 
  mWasabi fluorescence than with Coomassie stain---but this isn't an artifact.

  Perhaps using less mRNA gives the proteins more time to fold correctly, and 
  so ends up creating more functional protein.  I'd want to repeat this 
  experiment with more of a concentration gradient to verify this.  NEB does 
  say that the optimal amount of mRNA needs to be optimized for each gene, so 
  maybe that's just what I'm learning the hard way here.

- There's a slight discrepancy between the mWasabi fluorescence (A) and the 
  Coomassie stain (B).  With the former, the 1 µM mRNA template seems to 
  produce the most protein.  With the latter, the DNA template instead seems to 
  produce the most.  This may indicate that the SDS is affecting with mWasabi 
  fluorescence.  I don't know why SDS would have a different effect on the two 
  lanes, but I trust the Coomassie more.

- It's interesting that imaging Coomassie with NIR fluorescence (ex: 658 nm) 
  reveals bands that aren't visible colorimetrically (B).  This is pretty clear 
  evidence that NIR fluorescence is a more sensitive way to image Coomassie 
  staining, in agreement with [Butt2013]_.  This will definitely be good to 
  know going forward.

- mWasabi runs a little higher than it should, but I'm not going to read 
  anything into that.

Express mWasabi (gradient) --- 2020/02/26
-----------------------------------------
I'm going to do a serial dilution of mRNA from 10 µM to 0.1 µM in 7 steps.  I 
know that expression is higher at 1 µM than at 10 µM, so I'm hoping that 0.1 µM 
is low enough to see expression go back down again (since I want to find the 
maximum).   I chose 7 steps for several reasons:

- To get reasonably fine grained data.

- To include 1 µM, to compare to yesterday's experiment.

- To use exactly 2 aliquots of PURExpress (with 5 µL reactions).

.. protocol:: 20200226_serial_purex_page.txt

.. figure:: 20200226_express_mwasabi_mrna_gradient.svg

   (A) Direct imaging of mWasabi fluorescence.  The concentration refers to the 
   stock concentration of the mRNA template in step 2 of the above protocol.  
   (B) Coomassie staining.

- 1.0 µM mRNA template clearly optimizes protein expression.  It remains 
  surprising to me that protein expression *decreases* (rather than plateauing) 
  if too much mRNA is added, but the effect is very clear.

- In this case the mWasabi and Coomassie results agree.

- The Coomassie gel is warped because it dried out and started to roll up while 
  the image was being taken.

Express Zif268 (gradient) --- 2020/02/27
----------------------------------------
.. protocol:: 20200227_serial_purex_page_coom.txt

.. figure:: 20200227_express_zif268_mrna.svg

.. datatable:: 20200227_express_zif268_mrna.xlsx

   I couldn't identify a peak for Zif268, so I just measured the area of the 
   whole collection of peaks around 11 kDa.

- I can't see any Zif268 bands, but this isn't particularly surprising:
  
  - Zif268 overlaps with the dark bands at the bottom of the reaction, see 
    :expt:`35`.

  - mRNA templates may give lower expression than DNA templates, see the 2/25 
    result.  It's also possible that the Y-tag lowers expression even further, 
    since there's no stop codon.

- Gel densiometry doesn't reveal any clear trends in expression, although this 
  was a long shot anyways.  Don't read too much into the 10 µM lane being lower 
  than the others; due to a quirk in imageJ, this is the only one I had to draw 
  a manual background for.  I think that accounts for the difference.  It may 
  be significant that the 0 µM lane is lower than the others, though.

- I could use Ni-NTA to purify Zif268.  The construct I'm using (f11) has a 
  His-tag between the gene and the Y-tag.  

  This idea only makes sense because the bands that seem most responsible for 
  obscuring Zif268 do not to bind Ni-NTA (presumably indicating that they're 
  ribosomal proteins), see :expt:`34` and :expt:`33`.  One band near 11 kDa 
  remains after Ni-NTA purification, but it is relatively faint and seems to 
  run slightly below Zif268 (this is easiest to see in :expt:`34`).

  I could also use Strep-tag to purify Zif268, although I'd have to prepare a 
  new template (f25, f26).

  Either way, I'll want to start with DNA templates, just to make sure that the 
  purification works and I can see the protein.  Then I can do an mRNA 
  gradient.

- For the purposes of troubleshooting cDNA display, I could just keep using 
  mWasabi instead of trying to get Zif268 working.  I actually like this; it 
  helps focus on one problem at a time.  However, I would need to:

  - Clone mWasabi with Y-tag (and probably a His-tag, too).

  - Order and receive the new linker-N with Cy5.  The pseudo-linker doesn't 
    have puromycin, and the FITC-linker-N couldn't be used with mWasabi.

- I could've run the gel for 52 min without running Zif268 off the bottom.

Express Zif268 (StrepTag) --- 2020/02/28
----------------------------------------
.. protocol:: 20200228_purex_streptactin_page_coom.txt

.. figure:: 20200304_express_zif268_mrna_streptag.svg

- No protein was captured by the Streptactin purification.

- I did forget to add Zn to the protein expression reaction.  I wonder if that 
  could be the cause of the poor expression/purification.

- Maybe I should include a positive control next time.

Discussion
==========
- For mWasabi (f15), the optimal final concentration is 160 nM (which I 
  achieved by using 0.8 µL of a 1 µM stock in a 5 µL reaction).

- For Zif268, I was not able to detect protein expression.
