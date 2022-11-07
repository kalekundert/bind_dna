***************************
Check for denatured plasmid
***************************

[Hengen1996]_ suggests that plasmid preps can be contaminated by denatured 
plasmid.  This will manifest as a bright, fast-running band that can be 
eliminated by T5 exonuclease treatment.  I saw a fast band in one of my E-gels, 
but I've not captured any images of it as it's always run off the bottom.  I 
have T5 exonuclease, though, so I want to (i) capture the band and (ii) see if 
it's sensitive to T5 exonuclease.

.. protocol:: 20220921_t5_exo.pdf

.. figure:: 20220921_t5_exo.svg

- The T5 exonuclease does seem the be active (even though my tube was ≈1 year 
  expired).

  - The sr1 control gets noticeably fainter and more diffuse in the −enzyme 
     condition.

  - The 1124 plasmid has a very faint band around 500bp that disappears in the 
    −enzyme condition.  I think that is probably the denatured plasmid band 
    referred to by [Hengen1996]_.

- Denatured plasmid does not seem to be the reason for my low miniprep yields.

  - The pUC19 prep doesn't appear to have any denatured plasmid at all (i.e.  
    no bands are affected by T5 exonuclease treatment).

  - The 1124 plasmid has a very slight amount of denatured plasmid, but just 
    barely enough to be noticeable.  Plus, the 1124 plasmid seems clean by 
    Nanodrop anyways.
