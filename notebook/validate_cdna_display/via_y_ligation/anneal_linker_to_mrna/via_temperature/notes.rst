***************
Via temperature
***************

I was able to find a number of different temperature protocols for annealing 
oligos, see :expt:`15`.  They're similar for the most part---start hot, 
gradually cool---and I generally get the impression that any would be fine, but 
I want to test to see if there's any reason to prefer the longer, slower 
protocols.

Results
=======
Below are the controls I want to run:

- mRNA (f11) only, hold at 4°C
- linker-N only, hold at 4°C
- no PBS, hold at 4°C

Below are the temperature protocols I want to test (from :expt:`15`):

- Hold at 4°C.
- 95°C for 2 min, leave on bench until cool.
- 95°C for 2 min, 95°C→25°C over 45 min, hold at 4°C
- 95°C for 5 min, 95°C→25°C over 70 min, hold at 4°C
- 95°C for 5 min, 95°C→Tm at 1°C/min, Tm for 30 min, Tm→25°C at 1°C/min, hold 
  at 4°C

The Tm of the Y-tag sequence ``CAAGGGCGGGGGGCGGCGGGG`` is 88°C according to NEB 
TmCalc and 75–76°C according to SnapGene.  I'm inclined to trust SnapGene more, 
because I buy `their claim <https://www.snapgene.com/support/faq/>`_ that their 
Tm algorithm is more modern/accurate than NEB's, plus I'm not using Q5 here.

Based on :expt:`12`, I will use PBS for all reactions.

.. protocol::

   See binder, 2020/01/30.  By mistake I added 0.4 µL too much water to each 
   reaction (4.4 µL total), so everything is a bit too dilute.

.. figure:: 20200130_compare_gradient_annealing_10.svg

.. datatable:: compare_gradient_annealing.xlsx

   "High Band" and "Low Band" refer to the two highest bands in the above gel.  
   The saturated lower band is not considered.  The "intensity 10" image was 
   used for the analysis.

- I don't really know how good this assay is, because I don't have good 
  controls.  But it's at least a way to make a decision.

- I don't know why there are two upper bands, or which (if either) is the 
  intended "mRNA + linker-N" species.  Once I reorder linker-N with Cy5, I 
  could repeat this experiment and visualize the mRNA, which would probably 
  help determine what the two bands are.

- It's disconcerting that the 2 negative controls, "−PBS" and "Protocol A" (4 
  °C) performed better than most of the gradient annealing protocols.  This may 
  just mean that the assay isn't very meaningful, which wouldn't be all that 
  surprising.

- There does seem to be a correlation between longer incubation times and less 
  annealed product.  Maybe this is due to RNase contamination.  Again, if I 
  could see the mRNA, that might be informative.

- I decided to use "protocol B" (95°C for 2 min) moving forward.  I didn't feel 
  comfortable using either of the negative controls (since they're not supposed 
  to work), and B was the best of the rest.  It's also the most convenient, 
  which is nice.

