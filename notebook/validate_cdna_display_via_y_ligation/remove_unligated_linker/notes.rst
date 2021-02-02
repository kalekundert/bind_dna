***********************
Remove unligated linker
***********************
Most mRNA/cDNA display protocols include a step to remove unligated linker, 
although the specific step varies:

- [Keefe2001]_: TBE/urea PAGE
- [Cotten2011]_: LiCl precipitation
- [Barendt2013]_: MWCO spin column
- [Kubo2020]_: NaOAc precipitation

I can also think of some other approaches:

- No purification.  As long as the puromycin linker is not present in too great 
  of an excess, it might not cause problems.  It doesn't seem very reactive: 
  :expt:`52`.

- Silica column (e.g. RNeasy, Zymo)

- Drop dialysis

- Using o194 to denature unligated linker.  This could conceivably be used in 
  conjunction with any of the approaches above.

I can compare these different methods by PAGE, which will let me know:

- How well they eliminate the unligated linker.
- How much material they lose.
- How much the mRNA degrades.

Results
=======
2020/09/02:

This experiment compares using urea and o194 to denature unligated linker prior 
to doing a Amicon spin filter purification.

.. protocol:: 20200902_compare_o149_urea.txt

.. figure:: 20200902_compare_o194_urea.svg

.. datatable:: 20200902_compare_o194_urea.xlsx
   :sheet: densiometry

.. datatable:: 20200902_compare_o194_urea.xlsx
   :sheet: fold change

- The mRNA purified with o194 was considerably more degraded than the mRNA 
  purified without.  Possible reasons:

  - The heating step (95°C for 2 min) is damaging to the mRNA.  This seems 
    unlikely, though, since an identical incubation is used to anneal the mRNA 
    and the linker.

  - o194 is contaminated with RNases.  o194 is resuspended in EB, which is not 
    necessarily RNase free.  In fact, it could be relatively high in RNase 
    content just because I tend to use/handle that bottle when doing minipreps.  
    I could order more o194 and try resuspending it in RNase free water.  Or 
    maybe I could have IDT resuspend it; their solution is presumably 
    guaranteed to be nuclease free.

- o194 (7x) is more effective than urea (2x) at removing unligated linker, 
  although both are beneficial.  So it might be worth trying the above ideas to 
  see if I can use o194 without degrading the mRNA.

- The purification process seems to lose about half of the starting material 
  present in the annealing reaction.
  
  This assumes that I calculated the dilution for the annealing reaction 
  correctly and that I actually did recover 15 µL from the Amicon spin filter.  
  The latter I didn't explicitly check, but it's in the right ballpark: I used 
  4.4 µL for the gel and had about 10 µL left over.  The fact that the 
  intensities in the o129 and f85 lanes are about the same as the −ligate lane 
  also gives confidence that this is a meaningful comparison.

- The ligation yields observed in this reaction are consistent with previous 
  experiments.

- I don't know why the −filter lane didn't run cleanly.  Clearly the product is 
  fine, based on the lanes after it.  I've also run crude ligation reactions 
  before with no problem.  I must have done something wrong, although I'm not 
  sure what.

