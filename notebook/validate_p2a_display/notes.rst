********************
Validate P2A display
********************

Given that repA-display seems to be promising, I want to quickly see if 
P2A-display is equally promising.  The advantage of P2A display is that the 
coupling between the DNA and protein is covalent, which may allow more 
stringent treatments.  The downside is that P2A is larger, and might therefore 
exhaust the in vitro translation reaction more quickly.

Considerations
==============

P2A sequence
------------
The first step is to figure out the sequence of the P2A gene.  [Reierson2005]_ 
doesn't actually include this information, instead just stating that the y got 
the P2A vector from Isogenica Ltd.  This is definitely a red flag, but I'm 
still curious enough to try the idea at least once.  Interestingly, the 
Isogenica website now states that they make use of CIS display, so it's 
definitely possible that P2A display just doesn't work very well.

- I downloaded the complete P2 genome sequence from `NCBI Genbank`__.  I found 
  this entry via the accession number (AF063097) listed on the "Bacteriophage 
  P2" page on wikipedia.

__ https://www.ncbi.nlm.nih.gov/nuccore/AF063097.1

- Expression of P2A in lysate gave multiple product bands, possibly due to 
  internal start sites [Liu1993]_.  This is definitely something to be worried 
  about.

  - I ran the CDS through the Salis lab RBS calculator just the get an idea how 
    much internal expression there might be.

    .. datatable:: rbs_calc_max.xlsx

    - Both internal expression sites identified by [Liu1993]_ are predicted to 
      be more highly expressed than P2A itself.

    - There are also many other internal start sites even stronger than those 
      that [Liu1993]_ did not mention (although I think the sites mentioned by 
      [Liu1993]_ are backed up by in vitro expression data).

- A large C-terminal part of the sequence may be removable [Odegrip2001]_.

  On second thought, best to just use the whole protein.  [Odegrip2001]_ go on 
  to hypothesize that the C-terminal domain may be involved in cis-activity 
  (citing unpublished data).

- The F450/F454 mutant would be a good negative control [Odegrip2001]_.

Controls
--------
- With repA display, I used the following negative controls (at various 
  points):

  - STOP codon after mWasabi

    This was a bad idea, because I'm pretty sure it introduced an internal RBS.

  - Shuffling CIS/ori

    This was also a bad idea, because shuffling these sequences didn't seem to 
    fully eliminate their function.  Also, it wasn't totally clear to me where 
    either of these sequences were truly essential.

  - Deleting CIS/ori

  - Deleting repA

- For P2A display:
  
  - I can mutate P2A to eliminate the possibility of it covalently attaching 
    itself to the DNA.  Presumably it would still bind DNA, but I think it's 
    reasonable to expect that SDS PAGE would separate the two.

  - I can't eliminate the binding site, because it's located within the P2A 
    gene.

I think that the F450/F454 mutant should be the only negative control I need.

RBS
---
I'm wondering if it would make the most sense to use the native P2A RBS:

- The native RBS appears to be very weak.

- See :expt:`??` for a discussion of why a strong RBS might not be a good idea.
