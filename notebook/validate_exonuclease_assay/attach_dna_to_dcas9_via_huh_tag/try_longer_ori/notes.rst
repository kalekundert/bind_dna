**************
Try longer ori
**************

The primer I've been using to install the ori sequence (o102) has the minimal 
sequence needed to be recognized by PCV2-Rep.  There's plenty of evidence that 
this sequence should work (:expt:`46`), but Jorge also recommended adding 15 
arbitrary nucleotides between the ori sequence and the spacer and I want to 
give that a try.

Considerations
==============

Spacer sequence
---------------
Jorge didn't recommend a specific spacer sequence to use, he just said 
"NNNNNNNNNNNNNNN".  So I do need to pick a sequence.  Probably any choice is 
fine, but here are some options:

- Poly-A
- Random
- Random but controlling for GC content
- Corresponding sequence from PCV2 genome

  I found 4 PCV2 genomes after a very brief search.  All have the same ori 
  sequence and context.  I used the 12-nt ori sequence from [VegaRocha2007]_ 
  and included 6 and 15 nt of upstream context and downstream context, 
  respectively.

  - PCV2 refseq: cgctgtAAGTATTACCAGcgcacttcggcagcg
  - PCV2a:       cgctgtAAGTATTACCAGcgcacttcggcagcg
  - PCV2b:       cgctgtAAGTATTACCAGcgcacttcggcagcg
  - PCV2c:       cgctgtAAGTATTACCAGcgcacttcggcagcg

  After skimming some literature of PCV2 and looking at figure 7 from 
  [VegaRocha2007]_, I realized that it makes sense to include the palindromic 
  stem flanking the ori and the 2 hexamer repeats downstream of that stem.  
  Both features along with the ssDNA hairpin, are implicated in PCV2 binding:

  - gaagtgcgctgtaagtattaccagcgcacttcggcagcggcag

  It's also worth noting that the ori sequence used by Jorge is actually wrong: 
  the real ori begins with "tAAG" not "taAAG".

I'm going to try two sequences: 

- o204: the native genome context
- o205: exactly what Jorge told me (with N replaced by A).


Results
=======

.. figure:: 20201105_compare_ori_linkers.svg

- The highest MW bands (≈10kb) have both protein and DNA, but appear green 
  because the protein signal is very faint.

- The DNA is cleaner in this experiment than it has been previously, e.g.  
  :expt:`46` and :expt:`68`.  This is because I ordered new primers for this 
  experiment, specifically to get cleaner amplification.

- I don't know why the DNA with the A15 ori shifts dramatically in the presence 
  of dCas9-Rep without EDTA.

  - I've seen this same effect in previous experiments, both of which used ≈10x 
    less protein:

    - :expt:`46`: Nearly all DNA in shifted band.
    - :expt:`68`: Only about 20% of the total DNA is in the shifted band.

  The easiest conclusion to make is that the DNA is reacting with truncated 
  protein.  With the EDTA and −protein controls, it's hard to argue that Rep is 
  not involved in the shift.  And there is a *very* faint red band under the 
  shifted green band.

- I don't know what the ≈200 bp DNA species in the A15 lanes is.

  It's slightly bigger than f105, but I wouldn't rule out the possibility of 
  contamination.  This band didn't seem to be present in the E-gel I ran after 
  the PCR (not shown), so it probably isn't an incorrectly-amplified PCR 
  product.

