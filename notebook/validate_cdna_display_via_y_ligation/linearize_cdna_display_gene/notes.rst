***************************
Linearize cDNA display gene
***************************

The first step in cDNA display is to linearize the gene just after the Y-tag.  
This ensures that transcription stops such that linker-N can be attached to the 
mRNA by Y-ligation [Naimudden2016]_.  [Naimudden2016]_ linearized their DNA 
template by using PCR to add the promoter, the GS-linkers, the His-tag, and the 
Y-tag to the coding region.  I think that makes sense for directed evolution, 
where you'll only get out the coding sequence after each round.

In my case, i.e. transcribing from a plasmid and not doing directed evolution, 
it makes sense to include all of the aforementioned elements in the plasmid.  
This allows me to verify the sequence of all of those elements in advance, and 
prevents me from needing to use ridiculously long PCR primers.  The reverse 
primer [Naimudden2016]_ used must have been >100 nt long, and the forward 
primer wouldn't have been short either.

This leaves me with a few possible ways to linearize the gene used for cDNA 
display:

- PCR: This allows me to dilute the template, which may contain RNase A from 
  the miniprep.  The challenge is that the reverse primer would have to bind in 
  the Y-tag region, which is extremely G-rich.  The PCR product could be messy, 
  and I wouldn't want to use that.

- PCR with His-tag primer: The His-tag should be a good PCR primer site (50% 
  GC), will (probably) be present in all of my genes, and is near the Y-tag.  
  This primer would be 61 nt, which is on the long side, but should give better 
  amplification.

- XmnI digest: I designed the plasmid to contain a XmnI site just after the 
  Y-tag.  This would create high-concentration of homogeneously linearized DNA.  
  Unfortunately, XmnI would leave an extra nucleotide after the Y-tag.  I don't 
  expect that this would significantly affect Y-ligation, but I can't say for 
  sure.  This approach may also require additional purification of the DNA, to 
  remove RNase A from the miniprep.

- Type IIS digest: I could use a type IIS enzyme to create a cut right after 
  the Y-tag.  This would create a sticky-end, though, so I'd either have to add 
  blunting enzymes or read about how RNAP deals with ssDNA and position the 
  overhang such that the right product still gets transcribed.

Results
=======

XmnI digestion --- 2019/10/04
-----------------------------
.. protocol:: 20191004_prepare_dna_template_via_digest.txt

- The phenol/chloroform/isoamyl alcohol I ordered (Acros 327111000) is supposed 
  to be yellow.  I was worried when I noticed the color, because phenol can 
  turn yellow when it's oxidized, and oxidized phenol should not be used for 
  extractions because it can break the DNA backbone.  However, according to the 
  `product information page`__, the color is expected.  It is due to the 
  presence of 0.08-0.12% hydroxyquinoline, a "stabilizer" that helps prevent 
  the oxidation of phenol (described for a different product :download:`here 
  <product_info_phenol_equilibrated_stabilized.pdf>`).  

  __ https://www.fishersci.com/shop/products/phenol-chloroform-isoamyl-alcohol-25-24-1-stabilized-molecular-biology-dnas-acros-organics-3/ac327111000
  
.. figure:: 20191004_xmni_digest_49_51.svg

- I confirmed that all of these plasmids have only a single XmnI site, so I'm 
  really unsure why there are 3-4 bands for each construct, all of which seem 
  to be too big (although my EB ladder might be part of the problem).  Maybe I 
  should send the plasmid for NGS...

XmnI digestion --- 2019/10/11
-----------------------------
.. figure:: 20191011_xmni_digest_2_49.svg

- My plasmids really don't seem to be the size I think they are.  This has been 
  something I've observed consistently since starting here, but this result 
  brings it to the forefront of my attention again.

  There are only two explanations: My ladder is wrong, or my plasmids are 
  wrong.  Maybe the pUC19 I got with my MACH1 cells has a different sequence 
  than the plasmid I downloaded from SnapGene (or wherever I got that plasmid 
  map from).  If the sequence is wrong, it hasn't caused any problems so far, 
  but it'd be good to get that figured out.

  .. update:: 2019/10/24

     I sequenced the full p002 plasmid, and confirmed that it has the expected 
     sequence (i.e it is in fact 2686 bp).  

- This gels helps identify some of the bands from the 10/04 gel.  The bands at 
  3.0 kb seem to be the cleaved product, while the bands at ~2.5 kb seem to be 
  the uncleaved, supercoiled product (supercoiled because it runs faster than 
  the linear band).  I don't know what the bands at 5.5 and 3.5 kb are, 
  although one may be nicked plasmid.

  In any case, the 10/04 gel seems to indicate incomplete cleavage.  I thought 
  I'd calculated the appropriate amount of enzyme to add, but it seems that I 
  need more.  I should increase the enzyme to something like 2 µL, and pull 
  aliquots at timepoints until the plasmid is fully digested.

Y-Tag PCR --- 2019/10/11
------------------------
I thought it would be worth trying to linearize the cDNA display gene by PCR, 
because cleaning up plasmid DNA is difficult and I think I'd be able to test 
things faster using PCR.  It might also let me side-step whatever problem I'm 
having with XmnI digestion.

.. protocol:: 20191011_pcr.txt

   I used a gradient over 12 well (i.e. a horizontal gradient) and only used 
   the middle 8.  This gave me a more linear range of annealing temperatures 
   (see figure).

.. figure:: 20191011_amplify_ytag_primer.svg

- I only saw amplification with the longer primer, and only at quite low 
  annealing temperatures.  I wonder if the Ta predictions are less accurate for 
  such short primers.

- For the lanes with amplification, I see two bands.  One is the expected MW 
  (~410 bp), and the other is ~500 bp.  I don't know what the bigger band is, 
  but I'm not willing to use PCR unless it gives me clean product.

- I do wonder if it's possible that my PCR mix is going bad from being kept in 
  a not-very-cool refrigerator...

   - No, I've done plenty of successful PCRs since this.

- I might try designing a primer that anneals behind the Y-Tag.  I would still 
  have to use XmnI, but I wouldn't have to do phenol-chloroform extractions or 
  ethanol precipitations.  I would want the primer long enough that I can see 
  the difference after digesting it.  pUC-seq-ori (3) would work well for this: 
  it would give a ~500 bp product, which should go down to ~400 bp after 
  cleavage with XmnI.

Discussion
==========
- XmnI digestion seems to be an effective way to linearize the DNA.
  
- PCR primers binding in the Y-tag region do not amplify well.
  
- I haven't tested using a His-tag PCR primer or a Type IIS restriction enzyme 
  yet.  I'll revisit those ideas if there is evidence that the XmnI digested 
  template impairs Y-ligation.
