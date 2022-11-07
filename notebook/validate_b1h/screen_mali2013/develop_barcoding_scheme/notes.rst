************************
Develop barcoding scheme
************************

I need to find a way to (i) barcode the designed protein and the target 
sequence and ideally (ii) barcode a reference gene that I can use to normalize 
for cell population.

Approaches
==========

Nickase + polymerase
--------------------
.. figure:: schemes_nickase.svg

- This would require ≈10 bp barcode + ≈20 bp primer/restriction site in the DBD 
  oligo, which would take away from the protein sequence.

- I'd probably also need to add primers after the barcodes, to be able to PCR 
  the oligos.  It might be possible to 

- Could do 3-part assembly in the first step to add arbitrary sequence between 
  DBD and target site.

- Have to include either 5' or 3' end of DBD in the oligo.

Random barcode + hairpin
------------------------
.. figure:: schemes_hairpin.svg

- Can make a piece of DNA with a duplicated barcode by (i) ordering an oligo 
  with a hairpin and 20-30 nt or random sequence and (ii) adding polymerase.  
  It might also be smart to add a constant region after the random region so 
  that I can make linear dsDNA.

- I could ligate this double-barcode with the library oligos in a couple 
  different topologies.

- I'd have to sequence to figure out which barcode goes with which 
  protein/target pair, though.  And that might not be very tractable, because 
  I'd have to sequence the entire DBD.


Considerations
==============

Oligo vendors
-------------
- IDT oPool:
  
  - Up to 350 nt, no QC.  Similar to real oligo pools, but much smaller.
  - Order everything in one pool, then amplify to separate DBD from target.
  - 4-7 day turn-around
  - $109

- IDT Megamer: $600 for 500 nt: too expensive

- IDT Ultramer: Max length 200 nt: too short

- Agilent SurePrint (G7220A): Max length 230 nt, require at least 7.5K 
  sequences.

- Twist:

  - Up to 300 nt
  - No minimum size
  - 7-14 day turnaround
  - I'm kinda planning to use Twist long-term anyways, so maybe it would be 
    good to get started with them.

I like the idea of ordering from Twist, but I might also start with IDT 
oPools just for faster turnaround.

Arrangement
-----------
I can minimize the number of cloning steps I'll need to do by thinking about 
which components can/must be adjacent:

- Barcode downstream of Zif268

In order to minimize cloning steps, its important to put any parts that can be 
adjacent, adjacent.  The natural ways to do this:

- Barcode just downstream of the DBD gene.

- Barcode just downstream of the target site.

Deterministic vs random barcode
-------------------------------
[Boldridge2020]_ uses PCR to add a random barcode to each library member.  I've 
been avoiding this approach, because I thought I would just get too many 
different barcodes.  But I should probably not write it off.

Library members
---------------
- My ideal library:

  - ≈5 proteins and ≈5 target sites
  - Known affinities for every protein to every site.
  - A wide range of binding affinities.
  - At least 1 specific and 1 non-specific protein.

- [ElrodErickson1999]_

  - 5 proteins

    - All mutations between positions 18-24.

  - 3 target sites
  - Not every protein/target pair was measured.
  - No weak binders: all affinities in the range 0.1-70.0 nM
  - Method: EMSA

- [Yang1995]_

  - 3 proteins
  - 2 target sites
  - All 6 comparisons made
  - No weak binders: all affinities in the range 0.5-188.9 nM
  - Method: SPR

- [Mali2013]_

  - 16 proteins

    - Mutations span 63 aa
    - Proteins have 15-16 mutations each.
    - I'd probably limit mutations to one finger in my real designs.  But the 
      idea of using 60aa/180nt oligos seems like it could be good.

  - 16 target sites

  - All 256 comparisons made

  - The assay is very far removed from a direct measurement of binding 
    affinity, though: (i) express ZFs on cell surface, (ii) incubate with 
    labeled DNA, (iii) measure total fluorescent signal.

  - Additionally, these comparisons are very qualitative.  The underlying data 
    are not reported, so the only thing we have to go off is Fig 1d.  There are 
    only three distinct pixel values in that heatmap (dark blue, blue, and 
    light blue), so the most I can say is strong/medium/weak binding.

  - If I pick ZF01, ZF03, ZF09, and ZF10, I get a good mix of behaviors:

    - Strong specific binding
    - Medium specific binding
    - Non-specific binding

  - If I'm willing to measure binding affinities myself, then these library 
    members would meet all my "ideal" criteria listed above.
    
  - Regarding measuring binding affinities: I could probably do EMSA without 
    having to learn too much, or I could figure out how to do SPR/BLI...
