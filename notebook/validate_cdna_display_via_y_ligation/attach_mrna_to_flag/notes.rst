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

Results
=======

Optimize FLAG visualization
---------------------------
I realized that the FLAG peptide is small enough that it's kinda hard to see on 
a gel.  In this experiment, I tried two things to make it easier to see:

- Include vs. exclude the tRNA digestion step

- Bis-tris vs. tricine SDS PAGE

  - Tricine gels are known to give better resolution for peptides.

.. protocol:: 20210415_optimize_flag_visualization.txt

.. figure:: 20210415_optimize_flag_visualization.svg

Observations:

- If looking at the raw image for the tricine gels, the notched gel is the 
  10-20% gradient.

- I can't see the 17 kDa ladder band on the tricine gels.  I don't know why.

- The crystal violet loading dye moves into the gel, despite its positive 
  charge.  This is the case even with (i) only loading buffer and (ii) only 
  glycerol + crystal violet.  These lanes are not in the above figure, but are 
  visible in the raw data for the bis-tris gel.  I think this is probably due 
  to the SDS binding the dye.

  That said, crystal violet does seem to have less of a shadow that SERVA blue 
  G250/phenol red.  It is slightly fluorescent in the Coomassie channel, but 
  not at all fluorescent in the BODIPY/GreenLys channel.  It also doesn't seem 
  to quench green fluorescence like SERVA blue G250/phenol red does.

- FluoroTect GreenLys has a broad, low-MW smear.  I'm not sure what this 
  exactly is.
  
  - It's unaffected by RNase treatment.

  - It's somewhat fainter in the −RNase and +mRNA lanes.
  
  - It's narrower on the tricine gels than on the bis-tris gels.  The FLAG 
    peptide unfortunately falls within the smear in both kinds of gels, but 
    Zif268 would be easily outside the smear in the tricine gels.

- According to Thermo, the resolution ranges for the gels tested in this 
  experiment are roughly:

  - 16% tricine: 0-6 kDa
  - 10-20% tricine: 2-200 kDa
  - 4-12% bis-tris/MES: 3-250 kDa

- The +mRNA, −RNase condition has a faint band with slightly higher MW than the 
  tRNA.  I wonder if this is tRNA that's still attached to the peptide?

Conclusions:

- The Bolt gels actually look like the best option in this experiment.  But I'm 
  hesitant because they didn't look nearly so good in :expt:`99` (Apr 7, 2021).  
  Maybe that could be attributed to the loading dye, though.

- The RNase treatment doesn't help with visualizing FLAG.  It doesn't eliminate 
  the smear, and in fact makes it a little brighter, which makes the FLAG band 
  harder to see.  RNase treatment would be helpful if my product were just 
  slightly bigger, though.

- I don't think there's a clear reason to pick any of these gels over the 
  others.  I might revisit this once I start trying to attach mRNA.

