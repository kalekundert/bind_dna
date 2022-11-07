****************************
Debug NotI/HindIII digestion
****************************

While trying to clone p237 (:expt:`199`), I found that the NotI/HindIII 
digestion gave unexpected bands, and was not working well.  Here I attempt to 
understand why.

2022/10/21
==========
This is an analytical digestion, with the goal of understanding the bands in 
the gel from the 10/19 purification in :expt:`199`:

.. protocol:: 20221021_test_digest.txt

.. figure:: 20221021_digest_p236.svg

- The HindIII digestion is what causes the high-MW bands.

  - My only thought is that the bands are caused by HindIII not releasing from 
    the DNA.

    - I used loading buffer with SDS when I did the gel purification, and that 
      didn't affect these bands.  That doesn't really mean anything, but it is 
      the first way I would've tried to eliminate protein binding.

    - I could see what happens if I do a phenol-chloroform extraction, I guess.

    - If the bands are due to HindIII-binding, I probably don't need to worry.  
      All of the bands would be cleaved, so I could just purify the whole 
      thing.

  - NEB__ suggests that the appearance of unexpected high-MW bands in a gel may 
    be due to the enzyme binding the substrate (i.e. the same thing I 
    hypothesized above).  The recommended solutions are:

    - Lower the amount of enzyme used.
    - Add SDS to the loading buffer (which I already tried).

    __ https://international.neb.com/tools-and-resources/troubleshooting-guides/restriction-enzyme-troubleshooting-guide

- The uncleaved plasmid is pretty much indistinguishable from the cleaved 
  plasmid, so I basically just have to hope that the reaction goes to 
  completion.  Alternatively, I could clone in a larger drop-out region.

- It's a bit suspicious that the −enzyme control runs almost identically all 
  the +enzyme samples.  Circular plasmid, whether supercoiled or not, should 
  run differently than linear plasmid.

  - It might be worth doing a set of double-digestions with more obvious 
    expected patterns (e.g. 1 and 2 kb bands), to verify that both NotI and 
    HindIII are active.

  - XmnI should work well for this.

2022/10/25
==========
.. protocol:: 20221025_check_noti_hindiii.pdf 20221025_check_noti_hindiii.txt

.. figure:: 20221025_check_noti_hindiii.svg

  I forgot to run a ladder with this gel, so I copied a ladder from another 
  experiment (with the same gel and run parameters).  The bands should be 
  roughly in the right spots, but don't take them too seriously.

Observations:

- None of the enzymes (NotI, HindIII, or XmnI) are cutting efficiently.

  - The expected products are all present, except for the 67 bp one.  The PAGE 
    gel does a great job of resolving each of these bands.
  
  - Only a very small fraction of the DNA is cut twice, and none appears to be 
    cut three times.

  - Given that the XmnI digestion doesn't seem to go to completion, some of the 
    2.4 kb band may be the result of NotI/HindIII digestion without XmnI 
    digestion.

- Possible explanations for the poor activity:

  - My DNA doesn't have the right sequence.
    
    - Maybe there are just a lot of mutations, causing a lot of the restriction 
      sites to be wrong.  Or maybe I have a mixed population of plasmids.

    - I don't see any unexpected bands, which means that the issue probably 
      doesn't involve any big insertions/deletions, or contamination with the 
      wrong plasmid.

  - My enzymes are bad.

    - Expiration dates:

      - XmnI: October 2020
      - NotI: January 2024
      - HindIII: September 2023

    - Maybe age is a factor for XmnI, but it shouldn't be for the others.
    - Maybe there's something wrong with my −20°C?  Or the NEBNow freezer that 
      I got the NotI and HindIII from?

  - NEB troubleshooting tips for incomplete digestions: 
    https://international.neb.com/tools-and-resources/troubleshooting-guides/restriction-enzyme-troubleshooting-guide

    - "Slow sites": apparently some sites just don't work as well, and 
      incubation time needs to be increased to 1-2 h.

    - DNA contaminated with inhibitors: Doesn't say what these inhibitors are, 
      but does say that miniprepped DNA is susceptible to such inhibition.  Can 
      test for this by mixing a "control template" with the miniprepped DNA.  
      If there are inhibitors, the control won't be cleaved as efficiently as 
      it would on its own.

      - This idea does appeal to me, since I've been having so much trouble 
        with dirty minipreps.

- Can't calculate percent cleavage with this data, because the 2.4 kb bands (1x 
  cut) are very saturated.

Next steps:

- Same experiment, longer incubation.

- Same experiment, someone else's enzymes.

- Same experiment, with "control template".  What could that template be?

  - pUC19 from transformation kits?  Very dilute: 10 pg/µL.  Probably not 
    detectable.  Has XmnI and HindIII sites, but not NotI.

  - gBlock with AmpR?  PCR with cleanup?  Amplify with o2/o3, gel purify, use 
    that?  That would give nice fragment sizes (900, 700, 200, 70).  Wouldn't 
    be ideal for mixing with p236, but I could tell if the issue was due to 
    some sort of inhibition.

Discussion
==========
2022/10/28:

I was able to get colonies for p237, so I put this down for now.  I'll have to 
figure it out when I come back to do a real library, though.
