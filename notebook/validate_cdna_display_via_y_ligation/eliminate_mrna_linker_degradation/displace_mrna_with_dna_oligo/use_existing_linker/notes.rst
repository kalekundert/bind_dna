*******************
Use existing linker
*******************

See :expt:`-1`.  The simplest approach is to design an oligo to work with the 
existing linker, so I won't have to reorder any expensive reagents or reclone 
anything.

Considerations
==============

Oligo
-----
It's worth nothing that the RNA/DNA duplex could be broken up using either a 
DNA oligo complementary to the linker, or an RNA oligo complementary to the 
mRNA.  The advantage of using an RNA oligo is that I could use a longer 
toehold, and the concerns in the "Expression" section below wouldn't apply.  
The disadvantage is that I'd either have to order the oligo (very expensive, at 
least $217) or synthesize it (very short, hard to purify).

I'm going to start with a simple DNA oligo.  I'll come back to this idea if I 
need to optimize conditions.

Toehold
-------
I want to include a toehold so that the oligo binds the linker more tightly 
than the mRNA.  Specifically, this means that the oligo would include the 3 C's 
that make up the "arm" of the Y-ligation.

The specific oligo is o224.

Expression
----------
With the linker no longer annealing the mRNA, it's possible that region of the 
mRNA that use to anneal the linker will now be expressed.  This is a secondary 
concern, but maybe down the road, I could put stop codons in that region to see 
if they have any effect.

Results
=======

2021/02/08
----------
.. protocol:: 20210208_protect_mrna_with_o224.txt
.. figure:: 20210208_protect_mrna_with_o224.svg
.. datatable:: protect_mrna_with_o224_densiometry.xlsx

Observations:

- The ligation efficiency wasn't quite as good as it has been previously.  ≈27% 
  of the mRNA from this prep is attached to the linker, compared to ≈40% in 
  :expt:`1`.  This is still reasonable, though.

- The duplex oligo (o224) doesn't confer any protection against RNase H 
  treatment.  I can measure the fraction of the mRNA cleaved by RNase H using 
  both channels, and the results are pretty consistent: only ≈2% of the input 
  mRNA is not cleaved, regardless of how much o224 is present.  
  
- o224 does seem to affect how the cleaved linker runs.  Specifically, it seems 
  to make the linker run slower, which doesn't make sense to me.

  Did I order an oligo that anneals with the RNA by mistake?  No, I double 
  checked that o224 is complementary to o127 (a precursor to o129).

Conclusions:

- I think the most likely problem with o224 is just that it can't compete with 
  the mRNA in terms of kinetics.  In order for an free oligo to out-compete a 
  hairpin, the oligo would need much higher binding affinity (i.e. more 
  complementarity).  It'd also probably help to have a much more careful 
  cooling gradient, e.g. incubating for some time between the annealing 
  temperatures of the free oligo and the mRNA, and cooling slowly.

- I wonder if the fact that the linker has very low complexity (i.e. it's 
  mostly G) could explain why the Cy5 cleavage product runs slower when o224 is 
  added.  The idea is basically that o224 could bind off-register to the 
  linker, which might result in some mRNA remaining attached to the freed 
  linker.  If this were the case, it might help to replace the Y-tag with a 
  more reasonable sequence (e.g. an SR primer).

