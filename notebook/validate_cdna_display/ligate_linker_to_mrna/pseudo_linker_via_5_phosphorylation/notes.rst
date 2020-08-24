*****************************
Ligate via 5' phosphorylation
*****************************

The 5' end of linker-N must be phosphorylated in order for it to be ligated to 
the mRNA.  [Naimudden2016]_ used T4 kinase to phosphorylate linker-N during the 
ligation reaction.  Other Y-ligation methods, e.g.  [Kitamura2002]_, avoid the 
enzymatic phosphorylation step by ordering phosphorylated oligos.

I think that the latter approach (i.e. phosphorylated oligos) makes more 
sense.  This approach guarantees that every linker is phosphorylated, and 
eliminates the possibility of the phosphorylation step causing problems.  
The extra modification won't even add much to the price, given how expensive 
linker-N is already.  I will test this hypothesis by doing ligation 
reactions with posphorylated (o93) and non-phosphorylated linkers (o99), 
side-by-side.

Results
=======
.. protocol::

   See binder: 2020/1/31

.. figure:: 20200131_compare_phosphorylation_10_3.svg

   o93 has a 5' phosphate; o99 doesn't.  The " + ligation reaction" condition 
   includes T4 DNA ligase buffer, BSA, and incubation at 25°C and 65°C.  See 
   the protocol for more details.

- I don't know why the image quality is so bad.  Clearly something went wrong 
  with the gel imager.  I tried taking two images, both both had the same 
  problem.

- I wanted to quantify the progress of the reaction by gel densiometry, but I 
  think the poor image quality for make this analysis inaccurate.

- Both the phosphorylated and non-phosphorylated linkers appear to be ligated 
  very efficiently.  Looking at the intensity of the un-ligated linker
  bands, the reaction appears to have gone slightly further for the 
  phosphorylated linker.  This isn't surprising; I expected phosphorylation to 
  improve ligation.

- It's notable that the ligation efficiency in this reaction far exceeds the 
  efficiency I've seen so far with linker-N.  That could be because I used 
  optimized salt and annealing protocols for this reaction, but it may also be 
  that something about linker-N is more resistant to ligation.  I'll have to 
  repeat this reaction with linker-N to find out.

