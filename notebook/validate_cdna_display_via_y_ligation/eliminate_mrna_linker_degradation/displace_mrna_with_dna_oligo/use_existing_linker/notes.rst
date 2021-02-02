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
