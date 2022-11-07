***********************
Validate nickase scheme
***********************

One of the barcoding schemes outlined in :expt:`171` involves using nickase + 
polymerase to create the requisite barcodes.  My goal here is to find out how 
well this scheme works in practice.

Considerations
==============

RNAseq or qPCR?
---------------
- qPCR doesn't use barcoding, so I wouldn't really be testing that.  I'd be 
  testing whether or not (i) I clone the barcodes and (ii) the plasmid gives 
  good signal.  Both are valuable questions, but I think it would be better to 
  do RNAseq.

- I guess the drawbacks of RNAseq are that it would be more expensive and 
  slower.  But I need to get experience doing NGS, and I don't expect that 
  the MiSeq will be too expensive or too slow.

Samples to test?
----------------
- I have the 5 target sequences I've been using for my standard curve.  I'd 
  like to include some different Zif268 sequences as well, but maybe that can 
  wait.

- That said, I might argue that I'm not really testing the barcoding if I 
  don't have at least two proteins.  For example, if the protein barcodes 
  were somehow getting scrambled, I wouldn't notice unless I used multiple 
  proteins.

- [ElrodErickson1999]_ could be a good reference for this.  They have EMSA 
  data for 5 mutants and 4 target sites.  That's about the throughput I want, 
  and EMSA data is not the worst.

DBD orientation
---------------
- In order for the DBD and the reporter to be facing opposite directions, the 
  DBD barcode has to come after the DBD (i.e. in the 3' UTR).

- I guess that doesn't really matter, since the reporter barcode will be on a 
  short transcript, so will basically be at the 5' and 3' ends.

