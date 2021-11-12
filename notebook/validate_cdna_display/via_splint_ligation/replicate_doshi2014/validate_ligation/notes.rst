*****************
Validate ligation
*****************

Before doing the purmoycin coupling reaction from [Doshi2014], I want to make 
sure the splint-ligation protocol works.

.. protocol:: 20211109_validate_ligation.pdf 20211109_validate_ligation.txt

  I also tried running f144 on this gel, but the results were uninterpretable.  
  I don't think you can run dsDNA on a TBE/urea gel.

.. figure:: 20211110_validate_ligation.svg

- I think the ligation worked.  It's hard to be sure, because the 
  ligated/unligated bands only differ by 27 nt (the length of o278), but I 
  think there's a small shift between the +ligation and −ligation bands.

  My data also looks just like the data in [Doshi2014], Fig 3c.

- There are two bands in the purified product, but I don't think they represent 
  ligated/unligated product, because both are shifted slightly above the mRNA 
  band in the −ligase reaction.

- It would be easier to know if the reaction worked if the linker was Cy5 
  labeled.  Then I could just look for yellow bands.  It would also make it 
  easier to tell if the puromycin coupling reaction worked, for the next step.

- I just realized that the splint does not perfectly line up the 3' end of the 
  mRNA to the 5' end on the puromycin linker:

  .. literalinclude:: splint_overlap.txt

  Presumably this is intentional.  I'd have to look through the literature to 
  figure out the reason, though.

- The mRNA seems to have a ≈200 nt contamination.  I assume this is incomplete 
  transcript.  It seems to have been removed by the bead purification step, 
  though.

- Not sure if I should run the gel for more or less time:

  - I almost ran the splint off the gel, and I think I did run the linker off 
    the gel.  I'd probably see both bands if I only ran the gel for 30 min.  
    That said, both oligos are in 5x excess (IIRC), so I wouldn't really be 
    able to conclude that much about how far the reaction progressed from 
    seeing those bands.

  - The separation between the ligated/unligated bands is just barely 
    perceptible.  It might be easier to see the separation if I ran the gel for 
    70-80 min.

  If I had Cy5-labeled linker, I could run it for less and get all the 
  information I want.

- The linker-only control is not visible.  This may just be because I ran it 
  off the gel.  The splint-only control is faintly visible.  I'm surprised that 
  it's so much fainter than the reactions, though.  I tried to have roughly the 
  same amount in the controls as in the reactions.  I'll have to check my 
  calculations.

