********
PUREfrex
********

2021/04/21:

I'm coming back to this experiment after have done a lot of optimization with 
mWasabi, e.g. :expt:`101`.  The main takeaways from that work are that I 
should:

- Use PUREfrex rather than PURExpress.  I'm still not sure which version of 
  PUREfrex (e.g. 1.0 or 2.0) I should use, though.

- Use FluoroText GreenLys and TBE/urea PAGE to visualize Zif268 expression.

Considerations
==============

Which mRNA to use?
------------------
- RBS : Zif268 : His6 : Y-Tag

  - Already cloned (f11), so I can use it right away.
  - His6 is usable in PUREfrex (as opposed to PURExpress), but if I were 
    starting from scratch I might prefer an untagged protein.

- Y-tag vs GSGS

  - Y-tag:

      - Annealing sequence used by [Naimudden2016]_. 
      - 79% GC

  - GSGS:

      - Annealing sequence used by [Reyes2021]_.
      - 56% GC

  - I've used Y-tag in most of my experiments so far, including :expt:`101`.  

I'll start with f11.

mRNA concentration
------------------
See :expt:`99`.
