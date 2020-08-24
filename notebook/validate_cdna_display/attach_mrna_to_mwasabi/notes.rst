**********************
Attach mRNA to mWasabi
**********************
My goal for this experiment is to show that I can use cDNA display to 
covalently link protein to mRNA.  Note that this experiment will not include 
the reverse transcriptions steps of cDNA display.  After doing the expression 
reaction, I'll use SDS-PAGE to test if a covalent link was formed.  I'm using 
mWasabi because it will be easy to unambiguously visualize in the gel.  See 
:expt:`19` for a discussion of the expression protocol.

Considerations
==============

Controls
--------
- An ideal negative control for testing the puromycin coupling reaction would 
  be o129 without puromycin.  But that would cost $400, and would really only 
  be useful for this experiment.  I might not even need it (i.e. if the 
  reaction clearly works), so I'll proceed without for now, and think more 
  about ordering this oligo (o197) if I need it.

- In :expt:`19` (2020/02/18), my "controls" were to leave out reactions.  This 
  is nice because it also allows me to monitor the reaction and see where 
  problems arise (if there are problems).

  I do want to include a "−ligation, +expression" sample (which is not 
  represented in the aforementioned experiment, because I think it might be 
  informative to see mWasabi covalently coupled to just the puromycin linker, 
  but not the mRNA.

Volumes
-------
I'm trying to work out how much mRNA to load on the gel.  This is both to 
decide how much o194 to add the loading buffer, and the work out how big to 
make each reaction.  I'm using :expt:`19` as a reference:

- Assuming that I concentrated the mRNA to the maximum ability of the Amicon 
  spin filter (15 µL), I added 333 nM mRNA to the expression reaction.  Note 
  the ideal mRNA concentration as indicated by :expt:`18` is 1 µM, so I'm in 
  the right ballpark.

  .. math::

    333 nM = \frac{40 µL \times 125 nM}{15 µL}

- The concentration of the mRNA in the:

  - expression reaction: 26.6 nM
  - loading buffer: 12.4 nM

- The gel was overloaded, though.  In my notes, I thought it'd be better to add 
  only about half the sample, which would be 6.8 nM mRNA.

- I can do this reaction in just the same way, I just need to keep the volume 
  of the ligation reaction the same.

.. note::

   At this concentration, I don't know if I'd be able to see the mRNA even if I 
   stained for it.  But I'll keep the mRNA-only lane, because it could be 
   useful.

Ladders
-------
Since I'm not planning to stain my gel, it doesn't really make sense to run 
ladders (although I might anyways, in case I decide to post stain).  But for 
the future, I might also think about ordering a fluorescent protein standard, 
e.g. Invitrogen LC5928.  This is basically a normal ladder with AlexaFluor 
coupled to all the proteins.  It's kinda expensive ($400/125 µL), so I'll have 
to think if I'll get enough use out of it to be worth it.
