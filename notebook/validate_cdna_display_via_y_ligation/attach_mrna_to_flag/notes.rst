*******************
Attach mRNA to FLAG
*******************

2021/03/01:

Many mRNA/cDNA display methods are validated on the FLAG peptide.  This makes 
sense, since FLAG is probably very well-behaved by virtue of being short and 
very polar.  It should also be easy to detect by Western blot.

Given that FLAG seems to be a standard positive control and/or 
proof-of-principle, I think it would be smart to try getting to work in my 
hands.  My goal is to establish a good positive control, from which I can try 
to get display with Zif268 working.

Considerations
==============

Ordering
--------

- Display-F

  GTAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACCAATGGAC

- FLAG-control-R

  TTTTTCACCTGATCCGCTGCCTTTCTGTTTACCCTTGTCATCGTCGTCCTTGTAGTCCATTGGTATATCTCC 

- ???_REYES2021_LINKER

  /5Phos/CCCTTCACCTGATCCGCTGAAAAAAAAAAAAAAAAAA/iSp18//iSp18//iCy5//iSp18/CC/3Puro/

- ???_REYES2021_LIG_RT_16

  - IDT: "The internal version of this modification [iAzideN] is attached to 
    the oligo through a dT base. Incorporation of the internal version will add 
    a dT nucleotide at that position. To avoid adding an extra nucleotide, 
    replace an existing T nucleotide in your sequence with the required 
    modification."

  - The [Naimudden2016]_ linker has 6 nt 3' of the azide, which is consistent 
    with the use of 6 nt primers for RT.  Conveniently, this sequence has a T 
    in the right place to do the same.

  /5Phos/CCCTTCACCTGA/iAzideN/CCGCTG

- ???_REYES2021_LIG_RT_25

  - The [Naimudden2016]_ linker has 25 nt of complementarity (plus 3 nt of 
    overhang for Y-ligation).  It's possible that the extra length is 
    necessary, as the modification in the middle could interfere with base 
    pairing.

  /5Phos/CCCTTCACCTGATCCGC/iAzideN/GCCTTT
  21 nt

  /5Phos/CCCTTCACCTGATCCGCTGCC/iAzideN/TTCTGT
  25 nt

- ???_REYES2021_PURO

  o125 is almost identical to the puromycin arm of the linker used by 
  [Reyes2021]_.  The only difference is that o125 has an extra spacer-18.  I 
  think I can just reuse o125.

  o125:   /5DBCOTEG/AAAAAAAAAAAAAAAAAA/iSp18//iSp18//iSp18//iCy5//iSp18/CC/3Puro/
  reyes:            AAAAAAAAAAAAAAAAAA/iSp18/       /iSp18//iCy5//iSp18/CC/3Puro/


Conditions
----------
- F
  o

Protocols
---------
- Order everything

Preparation
===========
.. protocol:: 20210408_make_mrna_gel.txt

- I lost all the RNA on the spin column, because I forgot that the FLAG mRNA 
  only 40 kDa---not large enough to be retained by the 100 kDa MWCO filter.  

  For comparison, o237 is 17 kDa.  I don't think a spin filter will be adequate 
  to separate these.  For now I'll probably just skip this step, and maybe I'll 
  try to do a gel purification eventually.

