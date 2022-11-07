******************************
Validate single-barcode scheme
******************************

The idea behind this barcoding scheme is (i) to only have a single barcode 
identifying the DBD and the target and (ii) to normalize the barcode counts 
from the RNAseq data by isolating and sequencing just the plasmid DNA.  This 
differs from all of the other schemes I've considered so far, which focus on 
having matching barcodes appear twice in the plasmid, so that the normalization 
can occur entirely within the RNAseq data.

Advantages:

- This is very similar to the scheme I would use with the His+Ura constructs.  
  (The only difference is that there would be no need to express the barcode in 
  an auxotrophy selection.)

- The cloning should be significantly easier because there's no need to 
  manufacture duplicate barcodes.

- I can install constant sequences on both sides of the DBD.  I'm not sure if 
  all (or even any) of my other schemes can do that.

Disadvantages:

- The DNA and RNA samples have to be handled separately.  Although I can't 
  think of exactly how, this could potentially make the normalization less 
  accurate.  In other words, this isn't a true internal control, and internal 
  controls are very nice.

Notes
=====
- With this scheme, I don't have any choice w.r.t. where the barcode goes.  It 
  has to be downstream of the target-driven promoter.

- I can put the DBD on either side of the target, though.

- I might be able to save ≈40 nt on the DBD oligo by ligating a random barcode.  
  I won't worry about that for now, though, because it would add some steps and 
  it seems unlikely that 40 nt will make or break anything.

- [Boldridge2020]_ uses restriction sites to piece things together and PCR to 
  add the barcodes.  I'm 

Approaches
==========
- The ORI + DBD approach does make it easier to mess around with the target.  
  The target doesn't *require* any cloning steps, but I can add steps to (i) 
  install a reporter gene or (ii) shorten the oligo by installing the promoter.  
  This flexibility seems like a small but clear benefit.

ORI + DBD
---------
.. figure:: single_cloning_ori_with_dbd.svg

- I could clone a reporter gene (e.g. GFP) between the target promoter and the 
  barcode.  That might improve mRNA stability, and would make the target 
  fragment long enough to use Gibson assembly.  It would add another cloning 
  step, but it'd still be on a small library.

- Although it makes me uneasy, I've seen repeatedly that having the DBD 
  upstream of the AmpR gene, without a terminator, improves the assay.

ORI + target
------------
.. figure:: single_cloning_ori_with_target.svg

- Can use Gibson for the final step, because all the fragments are large.

