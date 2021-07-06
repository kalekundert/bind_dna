****************
Optimize inducer
****************

Lemo21(DE3) cells respond to increasing rhamnose concentration by decreasing 
protein expression.  This can be useful for tuning the expression of poorly 
soluble proteins.  

At this point, I don't know whether or not PCV2-Zif268 is poorly expressing.  
It is expressed in denaturing conditions (assuming its apparent MW is ≈8 kDa 
greater than expected): :expt:`113`.  This is my first time expressing it in 
native conditions.

Results
=======
2021/06/29:

.. protocol:: 20210629_optimize_inducer.pdf 20210629_optimize_inducer.txt

.. figure:: 20210701_optimize_inducer.svg

Observations:

- The 100 µM and 250 µM rhamnose conditions show a small amount of the ≈35 kDa 
  product that I observed in :expt:`113`.  This product is visible both in the 
  first elution lane and in the difference between the clarified lysate and 
  flow-through lanes.  I tried to quantify which condition had more of the 
  product, but the bands were too faint.

- The induced and non-induced cell pellet lanes were too gummy to load on the 
  gel.  I was hoping that pelleting a much smaller volume of cells than in 
  :expt:`113` would solve this problem, but it didn't.  Instead, I'll have to 
  either sonicate the cell pellets or treat them with DNase.

- I think I need to load more of each sample.  I adjusted to volumes in this 
  experiment to load a consistent amount of material in each well.  It's hard 
  to say if I achieved that goal, as there's no band that would be expected to 
  be present in each lane at constant intensity, but in any case the signal 
  across the board is too low.
  
- I don't like working with the agarose Ni-NTA beads.  They're hard to pellet.  
  I'm tempted to buy the Ni-NTA spin columns...

- I think the bright band at the bottom is lysozyme, which I added to help lyse 
  the cells.

Conclusions:

- PCV2-Zif268 does express solubly, albeit at low levels.

- A small concentration of rhamnose is beneficial for expression.

Next steps:

- Repeat with more sample, and maybe try more rhamnose concentrations

- Say it's good enough, do a larger scale purification, and try to get enough 
  for at least a few pilot experiments.
