***************************************
Measure Zif268 binding via [Noyes2008]_
***************************************

My goal here is to simply recapitulate the basic B1H assay using the controls 
from [Noyes2008]_, to make sure the strain and the plasmids work in my hands.  
To that end, I'm planning to individually transform and test target and 
non-target bait plasmids.

.. todo::

  - Make NM plates ([3-AT]: 0, 1, 2, 4 mM)
  - Make chemically competent cells
  - Transform p166+p167 (+ control), p166+p168 (− control)
    - Recover: LB + Kan + Amp
  - Pick colony, grow to ~OD 1
  - Wash
  - Plate dilutions


Preparation
===========
2021/07/15:

.. protocol:: 20210715_make.txt

No colonies.

2021/07/16:

.. protocol:: 20210716_debug_digest.pdf 20210716_debug_digest.txt

.. figure:: 20210716_check_ecori_noti.svg

Observations:

- Both enzymes appear to have full activity.

- The backbone I used for the ligation cloning reaction yesterday seems to be 
  fully linearized.  This makes me think that the ligation was the problem.

2021/07/16:

.. protocol:: 20210716_debug_ligase.pdf 20210716_debug_ligase.txt

.. figure:: 20210716_debug_ligase.svg

Observations:

- I tried to reuse an E-gel, and this is what I get.  It's not really clear to 
  me if the ligase is working well or not.  I'm going to repeat this experiment 
  more rigorously here: :expt:`124`

- Also, this wasn't a great experiment because the fragments didn't have to 
  reassemble into the original product.

2021/07/19:

I realized that I forgot to use PNK the first time I tried this reaction.  PNK 
is important, because the annealed oligos aren't phosphorylated.

.. protocol:: 20210719_make_p168.pdf 20210719_make.txt

Observation:

- I got colonies for the reaction with the intended amount of insert, but not 
  for the reaction with an accidental ≈20x excess of insert.  I'm a little 
  surprised this made a difference.  Maybe the problem is that not all of the 
  inserts were phosphorylated in the latter reaction.

2021/09/08:

- I tried to gel purify f131 and f132, but failed because the sample floated 
  out of the well.  I think I accidentally "eluted" in PE.  That said, the 
  Nanodrop concentration looked reasonable; not sure how to explain that.
