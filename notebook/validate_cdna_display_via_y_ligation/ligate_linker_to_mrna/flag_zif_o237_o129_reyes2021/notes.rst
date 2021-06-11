*********************************
FLAG/Zif, o237/o129, [Reyes2021]_
*********************************

In preparation for trying to work with FLAG mRNA to reproduce [Reyes2021]_, I 
want to compare the ligation protocol from that reference to my established 
ligation protocol.  I'm going to do the comparison for both the FLAG and Zif268 
mRNAs, which have substantially different linker annealing sequences, and 
therefore may perform differently.

Results
=======
.. protocol:: 20210512_ligate_reyes2021.pdf 20210512_ligate_reyes2021.txt

.. figure:: 20210512_ligate_reyes2021.svg

Observations:

- Both protocols give similar yield.

  This is a little surprising, because :expt:`12` made it seem like salt was 
  required.  But perhaps that's only true when using PNK.

- The FLAG mRNA is very faint, especially relative to the Zif268 mRNA (which is 
  about 3x longer).  Specifically, the FLAG signal is about 26x fainter than 
  the Zif268 signal.  I'm not sure what accounts for this difference.  Perhaps 
  the Zif268 sequence has more secondary structure.  I could try staining with 
  SYBR green II (I think I have some).

- Note that the linker-derived yield calculations account for the fact that the 
  two protocols have different mRNA:linker molar ratios.

- Not sure why there are 3 FLAG mRNA bands, or which is full-length.  The 
  control didn't show up; I'm guessing because it was just too dilute (in 
  conjunction with the overall faintness of the FLAG mRNA bands).

- Not sure what the faint GelGreen signal around 53 nt is.  It's present in the 
  ladder lane, which means it's in the loading buffer.  I suspect something to 
  do with o194, since the loading dye was crystal violet.

  .. update:: 2021/05/19

    The 53 nt band is due to the loading buffer (probably o194).  See 
    :expt:`115`.

Discussion
==========

- The various differences between then ligation protocol I came up with and the 
  [Reyes2021]_ ligation protocol don't seem to be significant.

- The advantage of my protocol is that the pipetting steps are easier/more 
  accurate, by virtue of including the ligase in a master mix.

- As I've seen before, the reaction only goes to about 50% completion.
