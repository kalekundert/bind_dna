**********************
Attach mRNA to mWasabi
**********************
My goal for this experiment is to show that I can use cDNA display to 
covalently link protein to mRNA.  Note that this experiment will not include 
the reverse transcriptions steps of cDNA display.  After doing the expression 
reaction, I'll use SDS-PAGE to test if a covalent link was formed.  I'm using 
mWasabi because it will be easy to unambiguously visualize in the gel.  See 
:expt:`19` for a discussion of the expression protocol.

.. toctree::
   :hidden:

   first_try/notes
   try_spacer18_linker/notes
   try_high_salt/notes

Brainstorming
=============
2020/09/09:

So far, I haven't seen any mRNA/protein coupling in any of my reactions.  But 
I'm convinced that this is possible, because there are many reports of 
successful mRNA display out there.   Some ideas for how to get the coupling 
reaction working:

- Why is my mRNA getting degraded?

  I've seen this in all of my reactions, and I really have no idea what's 
  causing it.  But whatever it is, it could be the reason for the coupling not 
  working (e.g. the mRNA gets chewed up before the puromycin can react).

  It's kinda weird that there's no smearing.  The Cy5 band looks exactly like 
  free linker.  If this is degradation, it's going all the way to completion 
  very quickly.  It's enough to make me wonder if the ligation actually failed, 
  although I have a lot of confidence in the assay showing that it worked.

  I've got a feeling that I'll need to figure this out before anything else.  
  Some ideas:

  - Use PUREfrex: See :expt:`67`

  - Timecourse, +/- all PURExpress components.

    - Stain with GelGreen to see if we can actually see the RNA, i.e. to 
      distinguish between the RNA being degraded and the linker being 
      not-ligated.  GFP will be in the same channel, but hopefully that won't 
      confuse too much.

- Use optimized C-terminal sequence [Nagumo2016]_.

  I'm not very enthusiastic about this:

  - The authors saw some coupling even in their initial construct, and then 
    just optimized it to about 40%.  In contrast, I don't see any coupling.

  - The optimal sequence could be scaffold-dependent, i.e. if I try it and it 
    doesn't work, I won't really learn anything.

  - The paper (like [Barendt2013]_) did not include a high-salt incubation.  
    Since the initial mRNA display papers so strongly emphasize the importance 
    of that step, I just get suspicious when it's missing.

- Overnight -20°C incubation

  This step is mentioned by a few protocols, including [Liu2000]_, [Ma2011]_.  
  I got the impression from [Liu2000]_ that the high-salt incubation made this 
  unnecessary, but I guess it's still worth trying.

- Rabbit reticulocyte lysate.

  This is the "standard" expression system for mRNA display.  Although I really 
  believe that PURExpress *should* work better, and *should* have less RNase 
  activity, it might be worth sticking more closely to the established 
  protocols.  
  
  To do this experiment, I'd have to clone the eukaryotic enhancer back into my 
  mRNA.  I'd also have to buy lysate, although it's pretty cheap (e.g. Promega 
  L4151).

- Ribosome is not stalling.

  This is hard to believe.  Ribosomes are known to stall at RNA/DNA junctions, 
  which my mRNA has.  It also may be stalled even earlier by the RNA/DNA 
  duplex.  On top of all that, I'm using the same sequences that presumably 
  worked for [Naimudden2016]_.  That said, [Nagumo2016]_ reported that rare 
  C-terminal codons improve coupling by helping stall the ribosome.
