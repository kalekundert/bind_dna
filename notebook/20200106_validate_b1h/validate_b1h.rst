************
Validate B1H
************

B1H is an established assay [Meng2005]_, [Noyes2008]_ that should be able to 
measure library-vs-library binding between designed proteins and DNA targets.  
Compared to the *in vitro* assays I'm also developing [:expt:`1 
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
  than the in vitro assays (e.g. factors like crowding on plates or differences 
  in plasmid copy number could be significant sources of variability).  But it 
  also may do a better job of amplifying the signal.  This could go either way.

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

  This might be less of a concern with the omega-based system [Noyes2008]_.

.. toctree::
   :glob:
   :hidden:

   /20200312_measure_zif268_binding_via_noyes2008/*
   /20200312_combine_b1h_plasmids/*

Considerations
==============

Nomenclature: "bait" vs. "prey"
-------------------------------
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

Selections: positive vs. negative
---------------------------------
[Meng2005]_:

  The HIS3 and URA3 reporter genes allow positive and negative selections to be 
  performed in a bacterial strain where the bacterial homologs are deleted.  
  Growth of cells on minimal medium containing 3-amino-triazole (3-AT), a 
  competitive inhibitor of HIS3, provides selection for an active promoter.  
  Growth of cells on medium containing 5-fluoro-orotic acid (5-FOA), which is 
  converted into a toxic compound by the uracil biosynthesis pathway, provides 
  selection against an active promoter.  Reporter vectors harboring a binding 
  site for the bait can be isolated by selecting for increased levels of HIS3 
  expression. Reporter vectors containing DNA sequences that activate the 
  promoter independent of the bait (self-activation) can be eliminated by 
  selection against URA3 expression. Thus, recognition sequences for the bait 
  can be isolated from the library of prey by a combination of positive 
  selection in the presence of the bait and negative selection in the absence 
  of the bait.

RNAP subunits: alpha vs. omega
------------------------------
RNAP is composed of 2 α-subunits (rpoA), 2 β-subunits (rpoB), and 1 ω-subunit 
(rpoZ).  Interestingly, the ω-subunit is not essential: it can be knocked out 
and cells will still grow normally.

In [Meng2005]_, the DNA-binding protein is fused to the α-subunit.  In 
[Noyes2008]_, the DNA-binding protein is instead fused to the ω-subunit, and 
the endogenous ω-subunit is knocked out.  The latter approach gives better 
sensitivity, presumably because there is no competing ω-subunit.  Both Scot 
Wolfe and Marcus Noyes recommend using the ω-based B1H system for most 
applications.  The only exception is for dimeric DNA-binding proteins.  These 
work better with the α-based system, presumably because RNAP has 2 α-subunits 
and only 1 ω-subunit.


Reagents
========

α-based B1H [Meng2005]_
-----------------------
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

ω-based B1H [Noyes2008]_
------------------------
- USO hisB- pyrF- rpoZ- (addgene #18049)

  - Tet resistance

- Reporter plasmid (pH3U3)

  - pH3U3-zif268 (addgene #18046)

  - This is very similar to the plasmid of the same name from [Meng2005]_.  The 
    only difference is a 14 residue deletion just after the Zif268 target site, 
    such that the target site is 10 rather than 24 bp upstream of the -35 box.  

    The spacing between the target site and the promoter is important, though.  
    See Fig. S6 of [Noyes2008]_.  The selection works the best when the target 
    sequence is 10/21 bp (ω-based system) or 14/25 bp (α-based system) upstream 
    of the -35 site.  These spacing requirements relate to the positioning of 
    the α/β/ω subunits in RNAP and the geometry of the DNA helix.  Small 
    changes in spacing are significant: the ideal spacing for the α-based 
    system appears to be ≈100x worse than ideal for the ω-based system.

- Bait plasmid (pB1H)

  - pB1H2w2-zif268 (addgene #18045): Zif268 with the lacUV5m promoter.
  - pB1H2wL-Prd (addgene #18040): TF (Paired) with the lpp-lacUV5 promoter.
  - pB1H2w5-Prd (addgene #18039): TF (Paired) with the lacUV5 promoter.
  - pB1H2w2-Prd (addgene #18038): TF (Paired) with the lacUV5m promoter.

  These plasmids are all based on the pB1H2 scaffold.  In [Meng2005]_, pB1H2 
  was only used in conjunction with pB1H1 (medium copy) for selections with 
  heterodimeric DNA-binding proteins.  In [Noyes2008]_, though, pB1H2 seems to 
  be the only plasmid used.  It probably doesn't make any great difference.  
  pB1H2 has a pBR322 origin (with the rop protein), which is similar in copy 
  number to pB1H1 (p15A origin).  Plus pB1H2 has AmpR, which is more 
  convenient.

  I'll need to order pB1H2w2-zif268 to use for my controls.  I thought about 
  also ordering pB1H2wL-Prd to get the lpp-lacUV5 promoter, but it will 
  probably be easier to clone myself if I want to do that.

Selective media
---------------
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


