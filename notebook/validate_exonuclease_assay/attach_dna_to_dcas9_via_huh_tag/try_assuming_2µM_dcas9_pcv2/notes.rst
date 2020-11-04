****************************
Try assuming 2 µM dCas9-PCV2
****************************

Based on the apparent excess of DNA in :expt:`68`  and the calculations in 
:expt:`70`, I'm going to try doing the reaction with the assumption that the 
dCas9-PCV2 fusion is not as concentrated as is seems.  Specifically, I'll 
assume it has an actual concentration of 1.84 µM based on the calculations in 
:expt:`68`.

Results
=======

2020/10/19:

.. protocol:: 20201019_2uM_cas9.txt

.. figure:: 20201019_2uM_cas9.svg

- There is a high-MW band in the +ori -EDTA reaction that matches the tagged 
  product seen in other gels, e.g. :expt:`46`, but it is a very minor product.

- For some reason the DNA channel didn't get imaged at all.  This was a fresh 
  batch of stain.  I tried a new staining procedure this time, so my guess is 
  that it didn't work.  Some ideas why:
  
  - The GelGreen stain was the last step, so it's not that the stain got washed 
    out by anything.
    
  - Note that the band representing the tagged protein should contain 
    covalently bound DNA, so the absence of signal can't be explained by the 
    DNA simply diffusing out of the gel.

  - SDS competes with DNA for GelGreen binding?  This doesn't seem likely.  The 
    protocol included washes to remove the SDS, and I've imaged DNA in SDS gels 
    even without washing before: see :expt:`46`.

  - The microwaving damaged the DNA?  Seems unlikely, DNA is not that fragile.
    
  I found a kit for doing EMSA experiments that comes with SYBR Green and SYPRO 
  Ruby for staining both DNA and protein: :download:`invitrogen_emsa_kit.pdf`.  
  The protocol for this kit calls for imaging the two channels separately, so 
  maybe there's not an easy way to get both channels at once.

- The protein is much more visible in this gel, while not being overloaded.  
  I'm still not confident in the actual concentration of the protein, though.



