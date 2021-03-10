***********************
Try SecM stall sequence
***********************

In ribosome display, a 17 aa sequence from the SecM gene [Nakatoga2002]_ can be 
used to stall translation and prevent the mRNA from dissociating from the 
ribosome [Contreras2007]_.  I haven't seen the same sequence used in mRNA 
display, but I thought it might be worth trying as a way of improving coupling 
efficiency.

I should note that using this sequence could also be detrimental.  
Specifically, it may prevent the mRNA from releasing from the ribosome even 
after the puromycin reacts.

Considerations
==============

Sequence
--------
- From [Nakatoga2002]_ abstract:

      Translation of SecM stalls unless its N-terminal part is "pulled" by the 
      protein export machinery. Here we show that the sequence motif 
      FXXXXWIXXXXGIRAGP that includes a specific arrest point (Pro) causes 
      elongation arrest within the ribosome.

- From [Contreras2007]_:

  - Peptide sequence: ``FSTPVWISQAQGIRAGP``

  - DNA sequence: not specified

- I used the primers referenced in [Contreras2007]_ to find the SecM peptide in 
  the MG1655 genome:

  - DNA sequence: 
    ``TTCAGCACGCCCGTCTGGATAAGCCAGGCGCAAGGCATCCGTGCTGGCCCTCAACGCCTCACCTAA``

  - Peptide sequence: ``FSTPVWISQAQGIRAGPQRLT*``

  - The [Contreras2007]_ sequence ends at the proline, but this sequence 
    matches perfectly up to that point.

- Final sequence, with TAA stop codon after proline:

  - DNA sequence: ``TTCAGCACGCCCGTCTGGATAAGCCAGGCGCAAGGCATCCGTGCTGGCCCTTAA``

  - Peptide sequence: ``FSTPVWISQAQGIRAGP*``

Stop codon
----------
[Contreras2007]_ claims to put a stop codon after the secM sequence, but I 
cannot figure out exactly how.  I tried to follow the cloning steps sepcified 
in the methods, but I did not end up with a stop codon.  [Nakatoga2002]_ claims 
to observe stalling with a TAA stop codon after the proline in question, so 
I'll probably just do that.

I thought about excluding the stop codon: the usual approach for mRNA display.  
But I'll keep it because it seems to get pretty stalled regardless, and I do 
want the protein to release after puromycin reacts.


