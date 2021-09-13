***************
Via refactoring
***************

My first approach to make a combined plasmid will be to aggressively eliminate 
any sequences that seem unnecessary.  This should leave me with the smallest 
possible plasmid, which should be easier to work with.

p169: optimize Ta
=================
2021/06/16:

I didn't get any colonies the first time I tried to make p169, so I tried 
optimizing Ta:

.. protocol:: 20210616_pcr.txt

.. figure:: 20210616_optimize_ta_p169.svg

Observations:

- The reaction amplifies well for a wide range of annealing temperatures.  I'm 
  going to continue using 59°C moving forward (the initial recommendation), but 
  I still don't know what the problem is.

2021/07/01:

I tried doing a large-scale PCR, to have clean product with which to try 
different ligation conditions.  However, the conditions from the above Ta 
optimization did not give clean product:

.. protocol:: 20210701_pcr_spin_cleanup_step.txt

.. figure:: 20210701_check_pcr_p169.svg

2021/07/02:

Based on the above, I repeated the Ta optimization:

.. protocol:: 20210702_pcr.txt

.. figure:: 20210702_optimize_ta_p169.svg

Observations:

- This time, none of the reactions gave the expected product.  Even the 
  relatively clean product in the 66.8°C condition is not the right size (≈2.5 
  kb instead of 5.8 kb).

2021/07/05:

Fitzy suggested that there might be something wrong with my water.  There might 
also be something wrong with my PCR master mix.  I'll setup four reactions to 
try all combinations of old/new water and old/new PCR mix.  I'll also use a 
fresh box of tips for everything.

.. protocol:: 20210705_pcr.txt

.. figure:: 20210705_compare_water_q5_p169.svg

Observations:

- The old PCR mix worked while the new mix didn't.  That doesn't really make 
  sense.

- The new water does seem to amplify better than the old water, although the 
  difference is slight.

- This time, the gel looks like what I saw on 7/1.

2021/07/05:

I repeated the temperature gradient with the new water and the old PCR mix:

.. figure:: 20210705_optimize_ta_p169.svg

Observations:

- This time, :math:`T_A = \pu{57°C}` seems optimal.  That seems plausible to 
  me, so I'll probably use that temperature going forward.

2021/07/06:

Given the stark difference between the old and new Q5 aliquots (and considering 
that both are form the same batch), I want to try using fresh Q5.

.. protocol:: 20210706_pcr.txt

.. figure:: 20210706_compare_q5_p169.svg

Observations:

- Neither master mix gave a significant amount of product this time.

2021/07/09:

I tried ordering new primers.

.. protocol:: 20210709_pcr.txt

.. figure:: 20210709_optimize_ta_p169.svg

Observations:

- Compared to the old primers, the new primers:

  - Give more of the intended product.
  - Require lower annealing temperatures.

- The product is still not clean.

.. protocol:: 20210714_pcr_gel_spin_x_gel_purify.txt

.. figure:: 20210715_gel_purify_p169.svg

Observations:

- The 4 and 2.5 kb minor products are the same ones I've seen previously.

- The PCR reaction worked very well this time, though.  The only thing I did 
  differently this time was to freshly dilute the template to 20 pg/µL from the 
  miniprepped stock.  I think it's very possible that my diluted stock was 
  somehow contaminated, and the cause of my problems.

- The gel purification yield was poor.  From a 50 µL PCR reaction, I recovered 
  10 µL of 22 ng/µL DNA.  That might even be an overestimate, because the 
  A260/A230 ratio was very poor.  Next time, I'll macerate the gel.

f83: optimize Ta
================
f83 didn't amplify well.  I initially tried to gel purify the correct band, but 
my yield was <10 ng/µL and the Golden Gate cloning failed.  Now I want to find 
PCR conditions that give better yield, so I can purify an appreciable amount of 
product.

2021/09/07:

.. figure:: 20210903_optimize_ta_f83.svg

- I ran the minor products off the gel, so I can't compare the purity of the 
  reactions.  Presumably, the higher annealing temperatures will have better 
  purity.

- I think I used :math:`T_A = \pu{60°C}` previously, so I can expect ≈2x more 
  product with :math:`T_A = \pu{67°C}`.

- I haven't seen the optimal :math:`T_A` yet, but this might be good enough.

2021/09/07:

.. protocol:: 20210907_make.txt

- Got undetectable yield.

- The ladder was much brighter than the sample (which is true in my :math:`T_A` 
  optimization as well).

2021/09/08:

.. protocol:: 20210908_pcr_gel.txt

.. figure:: 20210908_optimize_ta_f83.svg

- I ran the minor products off the gel again, for some of the lanes.  

- Annealing temperatures are inconsistent with the previous optimization 
  (2021/09/07).  This time, they're closer to the expected 60°C.  Presumably 
  this has something to do with the purity of the template.

- Yield seems much higher with the purified template, even though the 
  concentration of the template was too low to measure by nanodrop...

- I think it is important to gel-purify f81.

p170: Golden Gate
=================
2021/09/10:

.. protocol:: 20210909_make_p170.pdf 20210909_make_p170.txt

- Got 0 colonies.

  .. update:: 2021/09/13

    After the plate sat on my bench over the weekend, 4 colonies appeared.

- I don't think this plasmid is functional.  I can't think of any other reason 
  why the assembly would fail:

  - All of the fragments were purified and looked clean by nanodrop.
  - I added each fragment in exactly the recommended molar ratio.
  - I used reduced plates with less carbenicillin to account for the low copy 
    number.
  - I used a whole aliquot of fresh competent cells.

2021/09/13:

I did colony PCR on the colonies that appeared over the weekend to see if any 
are worth miniprepping and sequencing:

.. protocol:: 20210913_check_junctions.txt

.. figure:: 20210913_check_p170.svg

- None of the clones evince the expected bands.  All of the bands are about 2 
  kb too short, suggesting that the HIS/URA fragment was left out (somehow).

- The amplification is much cleaner for clone A, which is the only one of the 
  five that I miniprepped.  This is unsurprising, but going forward I think the 
  colony PCR results are good enough that I'd be able to detect a hit.
