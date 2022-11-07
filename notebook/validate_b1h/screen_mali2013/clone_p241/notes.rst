**********
Clone p241
**********

2022/10/05
==========

.. protocol:: 20221004_make_p241.pdf 20221004_make_p243.pdf 20221004_make_f195_f196.txt 20221004_make_p241.txt 20221004_make_f200_f201.txt 20221004_make_p243.txt
.. figure:: 20221005_electrotransform_p241_p243.svg
.. datatable:: 20221005_electrotransform_p241_p243.csv

- I got 10x more transformants this time than in previous transformations.

  - This was using 25 µL of cells that I'd already frozen and thawed twice, so 
    that means I don't need to worry about using larger volumes of cells, or 
    about using too many freeze/thaw cycles.

  - I think the reason for this is mostly that my DNA fragments were more 
    concentrated.

- My big plates were lawns.

  - Partly this is because I plated 2x more volume than previously

  - But mostly it's because my transformation was 10x more efficient.

  - The rule of thumb I found suggesting 10⁵ cells per 24 cm plate seems to be 
    pretty good.  Of course, it's not practical for big-big libraries, but 
    maybe I'll use 4 plates per library and hope it's good enough.

- I didn't need to dry my plates during the recovery.  Drying them in the 
  biosafety cabinet right after pouring them worked well.

2022/10/17
==========
I sequenced the "f202" fragment that I purified above, and the results suggest 
that things went really wrong:

- The fragment has the ori and Amp sequences repeated twice, but that same 
  sequence is on all of the plasmids used in these reactions, so there's noway 
  to say for sure where it came from.

  - It's an interesting idea to "barcode" each AmpR gene with som synonymous 
    mutations, to make it easier to track issues like this.  Probably more 
    effort than it's worth, though.

- The HIS, URA, and barcode sequences are all missing.

  - The lack of HIS/URA sequences means that f196 never got incorporated: 
    probably the last ligation reaction just re-ligated the backbone.

  - p240 is sequence-verified, so probably there was nothing wrong with f196.  
    I could sequence it to be sure, though.

- The fragment contains two BsmBI sites, suggesting that the digestion of f195 
  did not go to completion.  This may be because the site are too close 
  together, and the binding of one protein precludes the other.  Some ideas:

  - Re-order the library with more space between these sites.  I could do some 
    simple experiments to figure out the necessary spacing.

  - Do two digestions, with a protein-denaturing purification between them.  
    There are a couple of ways to do this, but probably they'll all have poor 
    yield.

- The junction sequence for the second BsmBI site has been changed from ATGG to 
  AGAT.  I can't figure out where the AGAT sequence came from.
  
  - It looks like I was thinking about using AGAT (found in sr132) in 
    :expt:`181`, but ultimately decided not to?

  - The sequence from the junction to the barcode—``AGATTGCCGGAG``— appears 
    just downstream of my "buffer" sequence in the plasmids that have it (e.g.  
    p182), although none of those plasmids were involved in these reactions.  
    p242 does use some of my "buffer" sequence, but it doesn't include this 
    part.  I also ordered it from a gBlock (I didn't clone it from one of the 
    plasmids that has this sequence), that isn't really an explanation.

  - The above sequence also appears in sr91 (also o105 and o166).  I don't 
    believe that any of those primers were used in any of these reactions, 
    though.

- My current hypothesis is:

  - p241 is either just f195 religated, or multiple copies of f195 ligated 
    together:

    - If p237 is not fully digested by BsmBI, it'll have a single sticky end 
      that could lead to polymerization.  (In contrast, a full digestion would 
      lead to incompatible sticky ends.)

    - NEB claims that BsmBI only needs a 1 bp overhang, and I should have much 
      more than that (the whole first site plus 2 extra base pairs).  But maybe 
      the first enzyme stays bound after making its cut, then gets in the way 
      of the second enzyme.

  - The fragments don't necessarily have two copies of AmpR and the ORI.  
    Instead, I purified a mixture of digestion products, and when the 
    sequencing program tried to align them, the sequences appear to be 
    duplicated.

    - After digesting "p241" with BsaI, I purified the ≈2.8 kb band.  The 
      sequencing result is 4.4 kb, though.  This means that the sequences must 
      not all overlap.

    - p237 is 2.4 kb and has three ≈evenly spaced BsaI sites, so the fact that 
      there was even a ≈2.8 kb band to purify means that the BsaI digestion did 
      not go to completion.  This is also supported by the fact that the f202 
      sequence has a bunch of BsaI sites.

- Next steps:

  - Need to use phosphatase on backbone to prevent re-ligation.  Should also do 
    −insert control transformations so I'll know if things worked.

  - Need to make sure digestions go to completion.

  - Decide whether to sequence p237, or remake it with phosphatase etc.

  - Remake p241 with phosphatase.

  - Sequence f196?

    - Two Sanger reactions to see both ends.

  - Do I need to sequence after each step?

    - I probably should for this screen, since I'm trying to debug things.  
      Once I have things working, I can probably stop.



