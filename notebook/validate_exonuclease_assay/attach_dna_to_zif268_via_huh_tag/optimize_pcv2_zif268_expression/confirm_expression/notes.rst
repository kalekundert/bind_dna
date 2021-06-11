******************
Confirm expression
******************

Before anything else, I want to make sure that my proteins are expressed.  
These experiments use small cultures and denaturing purification conditions, to 
maximize the chance that any expressed protein will be visible.

.. protocol:: 20210517_confirm_expression.pdf 20210517_confirm_expression.txt

.. figure:: 20210526_confirm_expression.svg

Observations:

- PCV2-dCas9 appears to be expressed, but at a relatively low level.  That's 
  ok, though; I can optimize expression later.

- The PCV2-Zif268 conditions have a bright ≈35 kDa band that looks like product 
  (e.g. present in the elution fractions), except that the product should be 28 
  kDa.

  The dCas9 elution fractions have a band at about the same MW, but it's much 
  more faint.

- The recommended loading volumes for the gel were not good.

  - The non-induced/induced controls were extremely gummy.  This may have been 
    a consequence of the way I froze and thawed the samples.

  - The flow-through, wash, and elution samples are very dim.  I should 
    probably add ≈5x more material next time.

  - I forgot to collect cleared lysate samples.

Conclusions:

- I think both proteins were expressed.

- I should stay vigilant about the PCV2-Zif268 construct, though.  Maybe think 
  about adding Zn²⁺ to the media?
