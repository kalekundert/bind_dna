*********************
Debug BsmBI digestion
*********************

In two cases (:expt:`202` and :expt:`204`), I suspect that I may be having 
problems due to the failure to fully digest two adjacent BsmBI sites.  Here I 
want to measure the progress of the digestion reaction in a more controlled 
environment.  Specifically, I want to use short fragments that can be readily 
visualized by PAGE.

2022/10/19
==========
.. protocol:: 20221019_check_bsmbi.pdf 20221019_make.txt 20221019_check_bsmbi.txt

.. figure:: 20221024_check_bsmbi.svg

.. datatable:: 20221024_check_bsmbi.xlsx

- The digest is inefficient, but the distance between the two restriction sites 
  doesn't seem to be an important factor.

- My starting templates aren't that clean.  It might be worth gel-purifying a 
  bunch of f206 to use in future experiments.

  - f206 because it starts off cleaner than the others, and to err of the size 
    of having too much space between the BsmBI sites.  Once I figure out what's 
    going on, I can try it again with f204.

2022/10/27
==========
I did a PCR annealing temperature gradient and a gel purification of f206.  
Including the results here in case I decide to make more later:

.. protocol:: 20221027_make_f206.pdf 20221027_make_f206.txt

.. figure:: 20221027_optimize_ta_f206.svg

.. datatable:: 20221027_optimize_ta_f206.xlsx

- :math:`T_A` does not affect the purity of the product.

  - I wonder if my template isn't pure.  Next time I run one of these gels, it 
    might be interesting to see what the template looks like.

- The reaction is pretty robust, but the empirical best :math:`T_A` looks to be 
  around 57°C.

- Total yield: 430 ng

  - Assuming that (efficient) PCR yields 50 ng/µL, and accounting for the fact 
    that 6 of my reactions weren't maximally efficient, I expected to get 4800 
    ng.  That implies 9% yield.

  - That said, I didn't measure my starting concentration, so I don't know for 
    sure what the yield was.

2022/10/31
==========
Try increasing the incubation time:

.. protocol:: 20221031_check_bsmbi_time.pdf 20221031_check_bsmbi_time.txt

.. figure:: 20221031_check_bsmbi_time.svg

.. datatable:: 20221031_check_bsmbi_time.xlsx

- I didn't add enough enzyme to the reaction.

  - The reaction is clearly working, in that I'm getting all the expected 
    products, it's just very slow.  

  - I didn't account for the fact that f206 has many more target sites per 
    nanogram than p237/p239.

    - 1 U of BsmBI is enough to digest 0.467 pmol of BsmBI sites in 1h.  (The 
      actual unit definition is based on λ DNA, which has 14 BsmBI sites and 
      MW=3.0×10⁷.)

    - Based on this, digesting 15 ng of f206 would require 1.27 U BsmBI.  
      However, I only used 0.15 U BsmBI in this reaction.

- I can't calculate how far the reaction progressed in the first two time 
  points, because the 82 bp band is saturated.  (The bands are not red in the 
  figure above because I did background subtraction, but they're saturated in 
  the original image.)

- I ran the gel slightly too long to see the 22 bp bands.  I'd probably be able 
  to see them if I ran the gel for 50m, but they don't really add anything to 
  the analysis.

2022/11/01
==========
Try using Esp3I instead of BsmBI, and also normalize the amount of enzyme 
relative to the number of target sites:

.. protocol:: 20221101_compare_bsmbi_esp3i.pdf 20221101_compare_bsmbi_esp3i.txt

.. figure:: 20221101_compare_bsmbi_esp3i.svg

- Both enzyme seems to work reasonably well given a 10x excess and a 1h 
  reaction time.

- It's odd that the 1x Esp3I reaction is significantly worse than the 1x BsmBI 
  reaction, but the 10x Esp3I reaction is slightly better than the 10x BsmBI 
  reaction.

  - Maybe I made a pipetting error when diluting the Esp3I?  I don't remember 
    noticing anything like that, though.

- Probably either enzyme is fine, I just need to be sure to to a long reaction 
  and a big excess of enzyme.
