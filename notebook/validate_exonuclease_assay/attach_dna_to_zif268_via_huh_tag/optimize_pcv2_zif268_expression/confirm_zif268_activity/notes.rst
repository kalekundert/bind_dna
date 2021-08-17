***********************
Confirm Zif268 activity
***********************

Before I can use my PCV2-Zif268 fusion to test the exonuclease binding assay, I 
need to confirm that the Zif268 domain is functional.

First try
=========
.. protocol:: 20210720_confirm_zif268_activity.pdf 20210720_confirm_zif268_activity.txt

.. figure:: 20210720_confirm_zif268_activity.svg

Observations:

- No protein is visible in the Coomassie channel.  This isn't surprising, 
  because I only added :math:`16 \times \pu{15 fmol} \times \frac{\pu{28000 
  fg}}{\pu{1 fmol}} \times \frac{\pu{1 ng}}{\pu{1e6 fg}} = \pu{6.72 ng}` 
  protein in the most concentrated lane.  For comparison, the detection limit 
  of Coomassie-stained gels is ≈10-100 ng.  Next time, I might start by picking 
  an amount of protein that should be easily visible, and scaling the amount of 
  DNA relative to that.

- A small amount of the DNA (≈11%) appears to be bound in the 16x and 8x excess 
  PCV2-Zif268 conditions.  This is noteworthy, because it means the protein has 
  at least some activity.

- There are two distinct "bound" bands, in the lanes that have any binding at 
  all.  Both are about equal in intensity.  I don't know the significance of 
  this.  In calculating the percentage of bound DNA, I included only the lower 
  of the two bands (not for any particularly good reason).

- I didn't load enough DNA.  I put a lot of thought into how much DNA I loaded, 
  and ultimately stayed pretty close to the 18.75 fmol recommended by 
  [Hellman2007]_.  But the signal is too faint.  I'd be inclined to load at 
  least 10x more DNA next time.  Note that I could get that by adding ≈3x more 
  DNA to the reaction and loading ≈3x more of the reaction on the gel.

Conclusions:

- I'm suspicious that my protein is not as concentrated as I think it is.  This 
  would explain why I see some activity at high "molar excesses" of protein, 
  but overall much less activity than expected.  I'm going to measure 
  concentration by PAGE/densiometry to check this.

  Alternatively, maybe my protein just isn't very active.  If this is the case, 
  the solution is probably to try optimizing the expression conditions.  I 
  might start by talking to Tina, since she told me that she has some 
  experience purifying Zif268.

Second try
==========
After observing that my PCV2-Zif268 was in fact much less concentrated that 
expected, I want to repeat the experiment with the intended quantity of 
protein.  See :expt:`121` for details on concentrating the protein.

.. protocol:: 20210816_confirm_zif268_activity.pdf 20210816_confirm_zif268_activity.txt

.. figure:: 20210816_confirm_zif268_activity.svg

Observations:

- The Zif268 domain is clearly functional.  

- Signal is too dim to quantify what fraction of the DNA is bound.
  
- The "bound DNA" band looks very different than it did previously.  In the 
  7/20 experiment, the bound bands were sharp and ran only ≈20% slower than the 
  unbound ones.  In this experiment, the bound bands were very diffuse and ran 
  ≈60% slower than the unbound bands.  (The unbound bands appear about the same 
  between the two experiments.)

  Looking back at :expt:`35`, the 7/20 experiment actually looks a lot like 
  free Zif268, while this experiment looks more like the fusions.  Or maybe 
  it's just a concentration thing, e.g. as the protein gets more dilute, DNA 
  that releases from a protein is less likely to be bound by another protein, 
  and so you get a tighter band.

- The 0.5x condition looks identical to all the other conditions, so I suspect 
  that the protein is at least 2x too concentrated.  
  
