*************
Try manganese
*************

[VegaRocha2007]_ shows that Mn is the best divalent cation for PCV2.  However, 
so far I've been using Mg because it's easier to find.  Now that I'm getting 
some unclear results, I want to try using the best cation to see if that 
improves anything.

Considerations
==============

Buffer composition
------------------
:expt:`46` has a list of several buffers that are relevant to this system.  
Some things to note:

- The actual buffering component is either Tris or HEPES in all of the above 
  buffers.  According the this `forum post 
  <http://www.protocol-online.org/biology-forums/posts/23999.html>`__, 
  phosphate buffers can sequester divalent cations and are not preferred for 
  enzymatic reactions.

- According to [VegaRocha2007]_: "There was no significant difference in 
  activity between pH 6 and 8, and comparable amounts of covalent adduct were 
  obtained at pH 6.4 and pH 8.4 (standard condition is pH 7.4)."

The specific buffers I'll be using are based entirely on what I can find in 
lab.  This ends up most closely following [VegaRocha2007]_:

Manganese buffer:
- 20 mM Tris, pH=7.0
- 100 mM NaCl
- 2.5 mM Mn(OAc)₂

Magnesium buffer:
- 20 mM Tris, pH=7.0
- 100 mM NaCl
- 2.5 mM Mg(OAc)₂

Negative control buffer:
- 20 mM Tris, pH=7.0
- 100 mM NaCl
- 2.5 mM EDTA

Results
=======

2020/10/15
----------

.. protocol:: 20201015_try_manganese.txt

.. figure:: 20201015_try_manganese.svg

  Note that this is a manually composited image from two separate scans, so the 
  channels don't line up perfectly.  The GelGreen channel in the second scan 
  was discarded, because it contained no signal.

.. datatable:: 20201015_try_manganese_coomassie.xlsx

  Densiometry results for the Coomassie channel.

.. datatable:: 20201015_try_manganese_gelgreen.xlsx

  Densiometry results for the GelGreen channel.

- The reaction seems to go to about 30% completion with both Mn²⁺ and Mg²⁺.

  - I don't have that much faith in the densiometry, because the bands are so 
    faint.

  - [VegaRocha2007]_ reports that the reaction goes almost to completion, so I 
    should expect better than this.

- The DNA appears to be in huge excess.

  - The "coupled" band (i.e. the highest one) contains about 30% of the protein 
    and about 3% of the DNA, suggesting that the DNA is in 10x excess.

  - Specifically, the excess is 11.94x, if you take the ratio of the means of 
    the two "% coupled" values for both channels.  If we assume that the DNA 
    concentration is correct, that would imply a dCas9-PCV2 concentration of 
    1.84 µM.

  - I was suspicious that the protein concentration was wrong, but it's not 
    clear whether or not this is the case: See :expt:`70`.  

- I don't know what all of the DNA bands represent.

  - Band 1 (>10 kb): DNA coupled to dCas9-PCV2

    - There's a matching band in the Coomassie channel

  - Band 2 (500-600 bp): PCR error

    - 500-600 bp in this experiment (a little imprecise because the ladder has 
      an unexpected band), where the intended amplicon is 272 bp.

    - 600-700 bp in :expt:`46` (imprecise for the same reason), where the 
      intended amplicon is 414 bp.

    - So about 200-300 bp longer than the intended amplicon.

    - I expect that one of the primers is just binding in the wrong place.  A 
      temperature gradient might help eliminate this band.

    - Strange that this band is missing from exactly one reaction.  I'm not 
      sure what to make of that.

  - Band 3 (400 bp): ???

    - Only present in conditions where the DNA should be coupled.

    - The Rep domain should cleave the ssDNA in the middle, but that band
      should be very small, e.g. ≈10 bp

    - [VegaRocha2007]_ does show that PCV2 Rep can catalyze the reverse 
      reaction, e.g. join two ssDNA segments with the appropriate sequence.  
      But the activity is inefficient, and I don't see how it could lead to 
      this band anyways.

  - Band 4 (300 bp): Uncoupled DNA:

    - It has exactly the expected size of 272 bp.
    - It also makes sense that f99 would be a bit shorter than f100, since it's 
      missing the ori sequence.

- The GelGreen seemed to be fully washed out by the Coomassie staining this 
  time.  Next time, I want to try the following protocol:

  - Wash the gel

    - 3x water + microwave

  - Stain with Coomassie

    - microwave + 5 min

  - Rinse Coomassie

    - 10 min

  - Stain with GelGreen

    - The GelGreen solution has salt, so this kinda mimics the second Coomassie 
      wash anyway.



Discussion
==========
I don't see a significant benefit to using the manganese buffer, but I'm going 
to keep using it for now.  It at least doesn't have BSA, which makes the gels 
cleaner.
