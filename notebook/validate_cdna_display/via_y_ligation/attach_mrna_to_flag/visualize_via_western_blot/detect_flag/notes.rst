***********
Detect FLAG
***********

After not seeing any binding in :expt:`126`, I need to troubleshoot my Western 
blot protocol before I can move on with optimization etc.  The specific things 
I want to try are:

- Only one condition: maximal antibody, incubation time, etc.
- Include dot blot, to rule out transfer problems.
- Load more FLAG, to rule out sensitivity problems.
- Larger antibody volume

It might not be possible to visualize the FLAG peptide.  [Reyes2021]_ doesn't 
show free FLAG being detected, only FLAG attached to the mRNA.  But I don't 
know if this is because (i) there was no free FLAG (seems unlikely), (ii) they 
just cropped the free FLAG band, or (iii) they couldn't detect free FLAG.

I might also think about using tricine gels in the future.  They should better 
resolve peptides.  For now I'm going to stick with Bolt, though, just because 
they're what I've been using and they're designed for western blotting.

FLAG peptide --- 2021/09/27
===========================
.. update:: 2021/09/30

   I realized today that I prepared the TBST incorrectly for this experiment.  
   I used 1 mL 10% tween instead of 1 mL 100% tween, so the detergent 
   concentration was 10x too low.  But I don't think this affects the results.  
   I got no signal, and if anything leaving out detergent should've increased 
   the background if anything.

.. protocol:: 20210927_debug_western.pdf 20210927_debug_western.txt

.. figure:: 20210928_debug_western.svg

- I see no signal for either the gel or the dot blot.  I already didn't really 
  think that the transfer was the problem, and this supports that.

- Reducing the wash time didn't solve the problem.  The background seems a 
  little higher than it did previously, though, and that could be due to the 
  less stringent washing.

- Using a greater volume of antibody didn't solve the problem either.

- The ladder lane didn't run very straight.  But this isn't the source of my 
  problems.

- I confirmed that the peptide I ordered (Bachem 4044101,
  https://shop.bachem.com/4044101.html) has the same 9-residue DYKDDDDK 
  sequence as expected by my primary antibody (Fujifilm 018-22381, 
  https://labchem-wako.fujifilm.com/us/product/detail/W01W0101-2238.html). 

I ordered another FLAG positive control, but I suspect the problem is either 
with the primary or secondary antibodies.  How do I troubleshoot that?

FLAG fusion control --- 2021/09/30
==================================
.. protocol:: 20210930_debug_western.pdf 20210930_debug_western.txt

.. figure:: 20211001_debug_western.svg

- The control lane has a clear band at the expected MW.  This indicates that 
  the primary and secondary antibodies are working fine, and that the problems 
  I've been having are specifically related to some aspect of the FLAG peptide 
  itself, e.g.:

  - Transferred through the membrane.
  - Not sticking to membrane.
  - Membrane binding interferes with antibody binding.

- I don't know why there are several bands in the control lane.  It shouldn't 
  be non-specific binding; the fusion protein should be the only thing in that 
  lane.

- I've noticed that the 8 kDa ladder band has been quite faint in most of my 
  western blots so far.  It was noticeably brighter in the experiment where I 
  only incubated the primary for 1h, but all the other ladder bands were 
  brighter in that experiment as well.

  I wonder if this is a symptom of running the transfer for too long.  Dima 
  recommended running the iBlot transfer for only 2 min, when working with a 1 
  kDa peptide.  The dot blot in the previous experiment was meant to test for 
  this possibility, but maybe it didn't work.  Dima was skeptical about it, and 
  thought that maybe SDS was necessary for FLAG to bind the membrane.

FLAG peptide + FLAG fusion control --- 2021/10/08
=================================================
Three changes:

- Tricine SDS PAGE (instead of bis-tris/MES)
- 2 min transfer (instead of 5 min)
- 1:10,000x primary antibody dilution (instead of 1:1000)

.. protocol:: 20211008_debug_western.pdf 20211008_debug_western.txt

.. figure:: 20211009_debug_western.svg

- There's no signal for the peptide.  With the positive control clearly 
  visible, this shows pretty convincingly that there's a problem with the 
  peptide.  I suspect that it's not adhering to the membrane.  
  
  It might be worth trying a nitrocellulose membrane, which may `bind better to 
  polar molecules`__.  Nylon membranes are also an option.  They're typically 
  used for nucleic acid blotting, but I think that that's because they have 
  amide groups that can crosslink with the nucleic acids.  I think they're less 
  polar than nitrocellulose.

  __ https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/western-blot-transfer-methods.html#blotting

- The 8 kDa ladder band is brighter in this experiment than it was in my 
  previous experiments.  I think this means that the shorter transfer is better 
  for peptides.  The high MW bands still tranferred very well though, which 
  means I probably still have room to decrease transfer time and/or voltage.

- The fusion control band is fainter in this experiment than it was previously.  
  This is probably because I used less primary antibody.  The "off-target" 
  bands is also fainter, but I still think that those bands are actually 
  tagged.

- Crystal violet has a strong fluorescent signal in the 700 nm channel (not 
  shown above, but clearly visible in the raw data).  It's not really a 
  problem, since I'm only using the 800 nm channel, but I think it'd be smart 
  to use a different loading buffer in the future.

Nitrocellulose --- 2021/10/18
=============================
.. protocol:: 20211018_debug_western.pdf 20211018_debug_western.txt


