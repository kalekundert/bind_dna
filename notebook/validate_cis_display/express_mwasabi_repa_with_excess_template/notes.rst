*****************************************
Express mWasabi-repA with excess template
*****************************************

2021/02/12:

For repA- and P2A-display, I only want to express one protein per template 
molecule.  This is very different from a normal translation reaction, where 
each DNA template would be used to make many mRNA transcripts, and each mRNA 
would be used to express many proteins.  Note that this amplification also 
helps slow down the reaction, which helps ensure that the reaction isn't 
exhausted before full-length protein can be expressed.

In order to get one protein per template, I either need to dilute the 
translation reaction or concentrate the template.  Diluting the translation 
reaction seems like a bad idea because it'll significant affect buffer and 
cofactor concentrations.  Concentrating the template is more plausible, but 
does present it's own issues.

Here I want to experiment with these ideas.

Considerations
==============
- I don't want multiple ribosomes expressing each transcript; since only one 
  protein can end up bound to the DNA.

- Maybe the best way to deal with that is just to titrate the amount of 
  template, e.g. with the goal of producing so much mRNA that only one ribosome 
  can express each one.

  - The final ribosome concentration in a PURExpress reaction is 2 µM 
    (:doc:`/appendix/appendix`).

  - PURE reactions have a 2.4:1 ribosome:RNAP ratio 
    (:doc:`/appendix/appendix`).  I couldn't find the same ratio for 
    PURExpress, but I assume it's similar.

  - NEB recommends a final concentration of 6 nM for the DHFR control template.  
    This is much lower than the final ribosome and (assumed) RNAP 
    concentrations.

  - The other important part of the equation is how much amino acids/ATP/etc.  
    are available for making protein:
    
    - NEB estimates that each ribosome is recycled 5 times with the control 
      DHFR template (:doc:`/appendix/appendix`).

    - P2A is about 10x longer than DHFR (≈2000 vs ≈200 aa), so I would expect 
      each ribosome to be used 0.5 times on average in these reactions. 

- I want to add a lot more template, but I also don't want to exhaust all the 
  amino acids before I have full length proteins.  I think that means I need a 
  much weaker RBS.

  - Note that I don't want a weak promoter, because having a small number of 
    transcripts would encourage each one to express multiple proteins.  What I 
    want is lots of transcripts, each expressed only once.
  
  - I'm also assuming that there's some kind of RNAP stalling mechanism to 
    enforce cis-activity.

- How can I increase the template concentration?

  - PCR is ≈50 ng/µL

  - The mWasabi-P2A gene would be ≈3.2 kb; MW≈1.93e6 Da

  - That works out to ≈25 nM.

  - I need to concentrate that ≈40x to get to 1 µM.

  - Options:

    - 50 µL PCR, ethanol precipitation, resuspend in 1 µL.

    - Plasmid prep, southern blot.


