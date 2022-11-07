******************
Develop qPCR assay
******************

Another way to which protein/target pairs are functional in a B1H assay is 
using RNAseq:

- Express a transcript with a barcode identifying the protein and the target
- Reverse-transcribe and sequence the transcript.

This does not require any selection, although it might be convenient with the 
URA3/HIS3 or GFP reporter.  It also promises high throughput, since I should be 
able to easily cover a library of ≈10⁸.

Considerations
==============

- How much do I want to adhere to the original plasmid?

  - Try to adhere to original plasmid as much as possible, and to be compatible 
    with all my insertion-site variants

  - Make new plasmid from scratch; use lessons from insertion-site variants, 
    but don't necessarily try to be over compatible.

  - It's good to have something that works, to iterate on.

GFP/RFP
-------
- I decided to merge this construct with the GFP/RFP construct, at least for 
  now.  It'll make cloning easier, and it might be useful to make side-by-side 
  comparisons of the two.

- Strictly speaking, I could just use the GFP/RFP cassette directly, and find 
  any necessary primers within those sequences.  But I decided to include 
  primer specifically for qPCR/RNAseq (SR001-4) to best reflect how I plan to 
  do this experiment for real.

Insertion sites
---------------
I currently have 6 single-plasmid constructs:

=========   ===========   =======  ========  ========
Construct   Zif268 Site   Barcode  Primer 1  Primer 2
=========   ===========   =======  ========  ========
p193        in f1         yes
p194        after f1      yes      SR013     SR071
p186        before AmpR   yes
p195        in f1         no
p196        after f1      no
p189        before AmpR   no
=========   ===========   =======  ========  ========

- I don't want to carry on with p193 or p195, because they didn't work well and 
  the way I designed this insert makes them unnecessary.

- Notes on primers:

  - p194:

    - SR013: This primer is between two opposing terminators installed with the 
      DBD insert.  The downstream terminator would be redundant with the RNAseq 
      insert I designed, so I want to remove it.

    - SR071: This primer is just upstream of AmpR.  It removes the rrnB 
      terminator protecting the reporter, which 

  - p186:

    - SR151: This primer is just upstream of the rpoZ-DBD lacUVm promoter.  It 
      removes a terminator upstream of that operon that I don't think is 
      important.  (If I wanted to keep the terminator, I could use SR123 
      instead.  Alternatively, I could add that terminator to the insert, to 
      help clamp down on expression.)

    - Something just downstream of the ORI.  I could use o345, which is 
      literally right after.  I could also use o284, which is what I used to 
      clone p194.

      I think I should go with o284 for now.  The point of this experiment is 
      to try RNAseq, not to minimize the plasmid as much as possible.

  - p196:

    - SR045: This primer is just downstream of rpoZ-Zif268.  The rpoZ-Zif268 
      operon doesn't have a terminator, but that's the whole point of this 
      construct.
      
    - SR071: See p194.

  - p189:

    - SR123: See p186, but this time SR123 (instead of SR151) is just upstream 
      of the rpoZ-Zif268 promoter, because the terminator isn't present.

    - Something downstream of the ORI: See p186.

Barcodes
--------
In order to have an internal control, I need two barcodes (one for the reporter 
and one for the internal control).  However, this may be very hard to clone:

- Throughput is limited by transformation efficiency, and each extra cloning 
  step I need will dramatically lower that ceiling.  If I need more than two 
  steps, I doubt that the assay will have sufficient throughput to be worth it.

- I need to install 4 sequence elements:

  - DBD coding sequence (≈300 bp)
  - target site (≈10 bp)
  - reporter barcode (≈20 bp)
  - internal control barcode (≈20 bp)

- I think I can do the assembly in two steps by including two promoters in 
  the oligo with the target site:

  .. figure:: cloning.svg

    Circle-within-circle: target site.  Arrow: promoter.  Straight line: 
    barcode. Blue block: DBD.

  Remaining cloning steps:

  - Linearize this circular construct between one barcode and the DBD.

  - Assemble into a plasmid.

  - Transform and isolate.

  - Linearize between the other barcode and the DBD.

  - Insert a fragment.

  - Transform and isolate.

  - Transform into assay strain.

- ≈250 bp would be necessary to encode all this: 10 bp target, 2x 40 bp 
  promoter, 2x 10 bp barcode, 50 bp terminator (not shown), 2x 20 bp PCR 
  primers (not shown), 2x 6 bp restriction sites, ≈40 bp insulation.  

  - The shortest [Chen2013]_ strong terminator is Bba_B0062-R (rrnC), which is  
    41 bp.  But I'm not clear if I should reverse it or not:

    - [Chen2013]_ gives it the "-R" prefixed, meaning it was tried in the 
      reverse orientation.

    - The sequence for this terminator in the supplemental materials table is 
      not reversed, though; is has the same sequence as BBa_B0062_ in the iGEM 
      database.

    - Bba_B0062 is itself the reverse complement of Bba_B0052_.  And I 
      confirmed by hand that Bba_B0052 is the natural orientation of the 
      promoter [Young1979]_.

    - BBa_B0052 is itself in the [Chen2013]_ data, and is in fact in the 
      "recombination-resistant medium-strength terminators" set.  It has a 
      measured strength of 16, compared to 110 for BBa_B0062.

    - Overall, I think this is probably a just bi-directional terminator.

  - The strongest [Chen2013]_ terminator is L3S2P21, which is 61 bp.

    - I'm already using this terminator elsewhere, so I'd have to refactor 
      first.

    - Not sure how precious space will be...

- The barcode would have to be in the 5' UTR.  I was hoping to avoid that, 
  since that makes it more likely that the barcode itself could affect 
  transcription, but it's unavoidable.  At least there will be a 20 bp PCR 
  primer between the promoter and the barcode.

.. _Bba_B0062: http://parts.igem.org/sequencing/part_analysis.cgi?part=BBa_B0062 
.. _Bba_B0052: http://parts.igem.org/sequencing/part_analysis.cgi?part=BBa_B0052

Transcription rate vs. concentration
------------------------------------
I want to measure the rate of mRNA transcription (:math:`\alpha`), since that's 
what's directly affected by DNA binding.  However, with qPCR/RNAseq I am 
actually measuring mRNA concentration (:math:`X`).  In order to relate these 
two quantities, I also have to know the rate of mRNA degradation 
(:math:`\delta`).  I don't know what the models for mRNA degradation are 
considered good, but I assume that it's reasonable for the rate to be 
proportional to the total concentration of mRNA (up to limits where the 
degradation machinery is saturated).  That gives the following model:

.. math::

  \frac{dX}{dt} = \alpha - X \delta

Assuming that steady-state is reached:

.. math::

  0 = \alpha - X \delta
  X = \alpha / \delta

This suggests that the mRNA concentration should be proportional to the rate of 
transcription, so long as the degradation rate is the same for all samples.  
Practically, though, there will be more mRNA if the degradation rate is lower, 
and that might give better signal:noise ratios.

mRNA degradation
----------------
According to [Selinger2003]_, mRNA in E. coli usually (but not always) degrades 
from the 5' end.  It might be interesting to do qPCR on both ends 
simultaneously, to see if I can measure a difference.

Protocols
---------
- I quantified RNA levels for my sgRNA project.  See the protocol here::

    ~sgrna/notebook/20180925_quantify_sgrna_levels/quantify_sgrna_levels.txt

- That worked very well, as I recall it, so I definitely want to start by just 
  doing the exact same thing again.

