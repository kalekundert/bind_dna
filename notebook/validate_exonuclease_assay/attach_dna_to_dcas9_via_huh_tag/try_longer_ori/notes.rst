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


