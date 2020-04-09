************
Validate B1H
************

B1H is an established assay [Meng2005]_ that should be able to measure 
library-vs-library binding between designed proteins and DNA targets.  Compared 
to the *in vitro* assays I'm also developing [:expt:`1 
<20190403_validate_cdna_display>`, :expt:`2 <20190419_validate_cis_display>`], 
this approach has some advantages and disadvantages:

Advantages:

- Established, and there for more likely to work.  That said, I will need to 
  make some changes in order to get a library-vs-library assay.

- More likely that the assay will support dimeric DNA-binding proteins.

- DNA binding proteins are not tethered to their targets.  This tethering has 
  the potential to affect the binding interactions.

Disadvantages:

- Lower throughput.  
  
  I'll be limited by how efficiently I can transform bacterial cells, and I'll 
  have to prepare my own competent cells because the assay requires a 
  non-standard strain (−hisB −pyrF).  My best competent cells preps have given 
  ≈1e7 transformants, so that's about what I expect the limit to be.

- Maybe less quantitative.

  In traditional B1H, DNA binding leads to expression, which confers a survival 
  advantage, which is measured.  This is a more indirect measure of binding 
  than the in vitro assays.  But it also may do a better job of amplifying the 
  signal.  This could go either way.

  I might be able to make this assay more direct by using a GFP or RNA-seq 
  based readout, rather than survival.  This would be a significant change from 
  the established protocol, however.

- Harder to control reaction conditions, e.g. salt concentration, competing 
  DNA, etc.  I'm worried that this might make it harder to tune the dynamic 
  range of the assay.

- Some DNA binding proteins may have toxicity: "The B1H system may not be 
  suited for determining the specificity of every DNA-binding domain. In 
  preliminary experiments, we observed that a bait containing the Max bHLH 
  domain was toxic in bacteria, and that selections using the bZip domain of 
  Giant did not yield a significant number of colonies." [Meng2005]_

.. toctree::
   :glob:
   :hidden:

   /20200312_combine_b1h_plasmids/*
   /20200312_measure_zif268_binding_via_meng2005/*

Considerations
==============

"Bait" vs. "Prey"
-----------------
The "bait" vs. "prey" nomenclature for B1H is really confusing to me.  The 
protein is the bait, and the DNA is the prey.  This is the opposite of how I 
think about it naturally, because I think of the protein as "going after" the 
DNA.  This is also opposite of the Y2H nomenclature, where the bait is the 
protein bound to DNA and the prey is the protein the recruits the polymerase.  

I think the reasoning behind this nomenclature is that you know what the 
protein is (so it's the bait) and you're trying to figure out what sequences it 
binds (so they're the prey attracted by the bait).  So the bait is the known, 
and the prey is the unknown.  Of course, this doesn't really apply to my 
application, where both are unknown, but that's what it is.

Reagents
--------
- E. coli hisB⁻ pyrF⁻ (addgene #12614)

   - Tet resistance

- Reporter plasmid (pH3U3)

   - pH3U3-mcs (addgene #12609)
   - pH3U3-zif268 (addgene #12610)

   - I might as well order the Zif268 one.  Both plasmids basically have the 
     same restriction sites, and the Zif268 one is a control I'll use.

   - Kanamycin resistance

- Bait plasmid (pB1H1)

   - pB1H1 --- rpoA fused to "dorsal RHR" (drosophilia) --- (addgene #12611)
   - pB1H1_zif268 --- rpoA fused to Zif268 --- (addgene #12612)
   - pB1H2 --- I think this is for heterodimers --- (addgene #12613)

   - I definitely want the Zif268 plasmid.  It has a different reverse 
     translation than the one I've been using, but I should be able to use it 
     right out of the box.

   - Chloramphenicol resistance

- NM and YM medium (His and 5-FOA selective, respectively):

   .. datatable:: media_reagents.xlsx

   .. note::

      I'm not sure why Met and Cys are excluded from YM.  I can't find 
      commercial media supplements lacking that specific combination of amino 
      acids, either.  −His, −Leu, −Trp seem to be the only common dropouts.

      I decided to just use a −His supplement, specifically "−His DO 
      Supplement" (Clontech 630415).  Especially given that yeast extract is 
      used in YM, it's hard to believe the experiment hinges on not having 
      Cys/Met.  And more amino acids probably just means faster growth.
      
      I chose this specific supplement because we already have it in the lab, 
      although I don't like that it's formulation seems to be proprietary.  For 
      example, these supplements often contain uracil and adenine, so I'm not 
      sure if I should add it exogenously.  I might think about switching to 
      `Sigma Y1751 
      <https://www.sigmaaldrich.com/catalog/product/sigma/y1751?lang=en&region=US>`_, 
      which does explicitly specify its composition.  I called Clontech to ask 
      about the composition and ended up giving them my contact info.  
      Hopefully they get back to me.

   .. note::

      The NM media recipe comes from [Joung2000]_.  Neither [Joung2000]_ nor 
      [Meng2005]_ seems to specify a carbon source...

      I think there is a standard recipe for "M9 minimal media" that includes 
      glucose as a carbon source.  From `Cold Spring Harbor 
      <http://cshprotocols.cshlp.org/content/2010/8/pdb.rec12295.short>`_:

      - 1x M9 salts
      - 0.4% glucose
      - 2 mM MgSO₄
      - 0.1 mM CaCl₂

      This is confusing because the specified MgSO₄ and CaCl₂ concentrations 
      differ from this standard recipe.  Perhaps that's why the salts were 
      specified but not the carbon source?  In any case, I think 0.4% glucose 
      is what I'd use.

   .. note::

      `This media <https://www.teknova.com/m9-minimal-medium-broth-csm-dropout-w-o-histidine.html>`_ 
      from Teknova is almost exactly what I want, except `CSM −His 
      <https://sunrisescience.com/shop/growth-media/amino-acid-supplement-mixtures/csm-formulations/csm-his-powder-100-grams/>`_ 
      (Sunrise Science, not Teknova) has 11 amino acids rather than 17, one of 
      which is methionine.

- Electrocompetent cell prep (not needed initially)


