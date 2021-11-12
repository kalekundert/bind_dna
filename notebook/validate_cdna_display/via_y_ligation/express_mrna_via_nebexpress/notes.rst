***************************
Express mRNA via NEBExpress
***************************

Before I an do experiments with NEBExpress, I need to determine the optimal 
amount of mRNA to use in reactions.  My plan is to closely follow the protocol 
from :expt:`18` (2020/02/26).

Results
=======

2021/02/17
----------
.. protocol:: 20210217_optimize_nebex.txt

.. figure:: 20210218_optimize_nebex.svg

.. datatable:: 20210218_optimize_nebex.xlsx

Observations:

- I see the highest expression at the highest mRNA concentration, but there are 
  diminishing returns.  In fact, the mRNA:expression relationship is pretty 
  linear up to 371 nM.

Conclusions:

- I think 400 nM is a good target concentration, because it should give 
  near-optimal expression without wasting mRNA.

- Unfortunately, the highest template concentration I can achieve with a 1 µM 
  stock is 240 nM.  I'd need a stock concentration of at least 1667 nM to reach 
  400 nM final.

Discussion
==========
- I didn't reach an expression maximum, but I think it's safe to say that 400 
  nM is a good template concentration.

