***********************
Validate hairpin scheme
***********************

The idea behind this barcoding scheme is to order the target and DBD oligos 
without barcodes, include a random barcode in the assembly, then sequence to 
map barcodes to target/DBD pairs.

Advantages:

- Most possible space allocated to library sequences in oligos.

  - This is the only way to fit Zif268 on a single 300 nt oligo.  Zif268 is 270 
    bp, and this approach requires only a primer on either end for PCR 
    amplification (15 + 15 = 30 nt).

Disadvantages:

- Need lots of sequencing to map barcodes to library members.

  - Note that this will not interfere with my ability to do small-scale 
    experiments.  The amount of sequencing needed is proportional to the number 
    of transformants, so I won't need more than 100,000 reads if I have ≈10 
    library members and keep ≈1000 transformants.

- So far I haven't been able to come up with a simple scheme for making the 
  barcode fragment.
