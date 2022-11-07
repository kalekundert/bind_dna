******************
Screen [Mali2013]_
******************

The goal of this experiment is to screen a small library of DNA-binding 
proteins/DNA targets using the growth assay (i.e. the His/Ura auxotrophy 
assay).  See :expt:`181` for a description of how the library was designed.  In 
addition to seeing if the screen can successfully identify known binders, I'm 
also hoping to get some experience with any pitfalls that might come up in 
cloning process.



Sequencing
----------
- f202:

  - Seems like f195 was not completely digested.  f202 still has both BsmBI 
    sites.

Notes
=====

Bioassay dishes
---------------
"BioAssay dish" seems to be the proper name for the large dishes used for 
plating colonies.  I found two manufacturers:

- `Corning™ Untreated 245mm Square BioAssay Dishes <Corning™ Untreated 245mm 
  Square BioAssay Dishes>`_
- `Nunc™ Square BioAssay Dishes 
  <https://www.thermofisher.com/order/catalog/product/166508>`_

Both of these plates are polystyrene, which means that they cannot be 
autoclaved and are therefore disposable.  I was not able to find anyone selling 
reusable plates (e.g. made of glass or more robust plastic).

Both plates are about the same price in general, but the Corning ones seem to 
be significantly discounted in B2P, so I ordered the Corning ones.

Gel purification
----------------
See :expt:`190` for more info.

Cloning improvements
--------------------
- Use Tina's plasmid instead of pUC19.  I expected pUC19 to give me good 
  miniprep yields, but it just doesn't.

- Add XmnI sites on either end of the barcodes.  That way, I can do NGS without 
  having to amplify the plasmid.

- Use phosphatase on backbone.

- Use B1H backbone for the target plasmid (p241) instead of the DBD plasmid 
  (p243).  This will allow me to remove some auto-activating sequences from the 
  library.

- Practically difficult to gel purify some fragments, e.g. f190.  Might be more 
  reliable to use PCR when possible, e.g. the early cloning steps.
