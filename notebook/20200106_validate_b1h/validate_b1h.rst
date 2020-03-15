************
Validate B1H
************

B1H is an established assay [Meng2005]_ that should be able to measure 
library-vs-library binding between designed proteins and DNA targets.  Compared 
to the *in vitro* assays I'm also developing [:expt:`1 
<20190403_validate_cdna_display>`, :expt:`2 <20190419_validate_cis_display>`], 
this approach has some advantages and disadvantages:

Advantages:

- Established, should be easy to get working

- Can use dimeric DNA-binding proteins.

- DNA binding proteins are not tethered to their targets.  This tethering has 
  the potential to affect the binding interactions.

Disadvantages:

- Lower throughput.  I'll need to transform two libraries.  I can try either 
  simultaneously transforming two plasmids, or remaking the competent cells 
  after the first transformation.  Either way is likely to give worse 
  efficiency than companies typically advertise, probably no better than 1e7 
  transformants.

- Maybe less quantitative.

  In traditional B1H, DNA binding leads to expression, which confers a survival 
  advantage, which is measured.  This is a more indirect measure of binding 
  than the in vitro assays.  But it also may do a better job of amplifying the 
  signal.  Care calibration will be required.

  I could possibly make this a little less indirect by using a GFP or RNA-seq 
  based readout, rather than survival.  This would be a next step, though.

- Harder to control reaction conditions, e.g. salt concentration, competing 
  DNA, etc.  

- Some DNA binding proteins may have toxicity: "The B1H system may not be 
  suited for determining the specificity of every DNA-binding domain. In 
  preliminary experiments, we observed that a bait containing the Max bHLH 
  domain was toxic in bacteria, and that selections using the bZip domain of 
  Giant did not yield a significant number of colonies." [Meng2005]_

In this experiment, I will measure the binding of my Zif268 controls using the 
B1H protocol as described by [Meng2005]_.

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

Consolidate plasmids
--------------------
[Meng2005]_ uses separate plasmids for the protein "bait" and DNA target 
"prey".  However, this approach doesn't work if both components are libraries, 
because there will be no way to tell which protein went with which DNA target.  
For my assay, then, I will need to consolidate both targets onto a single 
plasmids.

One way to do this is to join the two plasmids in vivo, by expressing a 
recombinase and including recombination sites such that the two barcodes end up 
next to each other in the recombined product.  Gleb and Pierce presented this 
idea in lab meeting about a year ago, so I should talk to them if I want to 
seriously pursue this.  Frankly, though, it didn't seem very robust.

The other way to do this is to clone everything into a single plasmid.  This 
also introduces some complexities:

- The size of the plasmid.  Larger plasmids are a little harder to work with 
  and are not transformed as efficiently, but I don't think this will be a real 
  problem.  The consolidated plasmid would be about 6 kb:

  - >2 kb: HIS + URA
  - >1 kb: rpoA + Zif
  - >1 kb: AmpR
  - >1 kb: ORI

  This would be bigger than most of my plasmids, but not really outside the 
  realm of the ordinary.  Plus, it'll probably be more efficient to transform 
  one big plasmid than two little ones.

- Gene direction.

  Initially I thought I would have the genes pointing away from each other, so 
  that run-on transcription wouldn't be an issue.  But this turned out the be a 
  bad idea, because it put the two promoters right next to the target, so they 
  would presumably both be affected by B1H.  (Orienting both genes in the 
  opposite direction so their promoters are far apart is not an option, because 
  then there's no tractable way to clone protein/target barcodes.)

  So instead, I'm orienting the genes back-to-back, and just separating them 
  with the strongest terminator from [Chen2013]_.  I might also want to include 
  another terminator just to really stop things (e.g. the strongest natural 
  terminator, to minimize sequence homology.)  The fact that I have positive 
  and negative selections gives me comfort that I'll at least know if the 
  terminator isn't working well enough.

- f1 origin

  pH3U3 has an f1 origin, which allows the plasmid to be packaged into phage.  
  Since this is not part of the B1H protocol, this sequence can be removed.

- Target is in MCS placed upstream of gene.

  I'm gonna replace that with PCR primers/Golden Gate junctions.

  - THe MCS doesn't work anymore anyways, most of the sites are no longer 
    unique cutters.

  - It really isn't that convenient to have a MCS.  PCR and GG are both easier.

  - The GC content is really high as it is, which would make PCR hard.

  - Because the authors just put the MCS there, it's unlikely that there's 
    anything really special about the sequence that needs to stay as it is.  
    I'll keep the spacing and hope that nothing else was important.

    Hopefully the very G-rich MCS isn't contributing to the poor promoter 
    strength, by making the duplex hard to open or something.  

- Include rrnB T1 + T2 before Zif268 target to insulate AmpR expression from 
  B1H expression.  Also insulates HIS/URA from cryptic promoters, although this 
  is unlikely to be a real problem.

- Cleanup:

   - Remove BbsI and BsaI sites.
   - 


- The plasmid copy number.  Currently the protein plasmid is medium-copy (p15A) 
  and the DNA target plasmid is low-copy (pSC101).  I can imagine that getting 
  the assay to work well may require tuning the copy-number of the plasmid (via 
  the origin) and the expression levels of the DNA-binding protein and the 
  reporter genes (via their respective promoters).
  
  To make this easy, I should design the plasmid to support modular Golden Gate 
  assembly.  It might even be worth striving for partial compatibility with as 
  established system, e.g. MoClo.

  .. note::

      I just brushed up on MoClo.  The original MoClo system [Weber2011]_ is 
      designed for eukaryotic genes, but two E. coli part libraries have been 
      described and made available on AddGene.  The first is CIDAR 
      [Iverson2016]_ and the second is EcoFlex [Moore2016]_.

      CIDAR mostly use the same overhangs as the original MoClo, but not with 
      the same meanings.  (So MoClo and CIDAR parts are not compatible, but 
      that's fine, they're meant for different organisms anyways.)  It's not 
      clear to me how CIDAR transcriptional units (TUs) are assembled, but 
      presumably I'm just missing something.  AddGene has both a CIDAR kit and 
      a CIDAR extension kit, which total to more parts than EcoFlex has.

      EcoFlex uses completely different overhangs than MoClo.  It also has 
      support for N-terminal tags.

      I don't think I can directly use CIDAR/EcoFlex, because I want my genes 
      pointing in opposite directions.


  

- Might be worth sending an email to the authors about this, just to ask is 
  they have any words of caution.

To address this, 



That works when only the bait is a library, but it my case I want both the bait 
and prey to be libIn my case, where both 

- Different copy

   - MoClo to expt with ORI

- Switch to Carb, fuck this Kan/Chlor shit

- Point genes different directions.  Don't want to get feedback loops

- Size
   
   Probably 6 kb total.  Not too bad.

In order to know which protein bound which target, I need to have the protein 
and DNA barcodes end up on the same DNA molecule.  Most likely 

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

Next Steps
==========
- Find:

   - Kan
   - Chlor
   - Unselective plates

- Prepare competent cells, maybe with bait already in there.

.. .. toctree::
      :glob:
      :hidden:

      /20190430_create_minimal_cloning_vector/*
      /20190603_express_zif268_in_vitro/*
      /20190614_optimize_rbs/*
      /20190625_purify_zif268_repa_via_reverse_his/*
      /20190626_purify_zif268_repa_via_ribosome_pull_down/*
      /20190626_purify_zif268_repa_via_rrna_digestion/*
      /20190711_purify_zif268_repa_via_affinity_tags/*
      /20190828_purify_zif268_via_imac/*
      /20190627_confirm_cis_display_via_labeled_dna/*
      /20190723_confirm_cis_display_with_fluorescent_protein/*
      /20190723_confirm_zif268_emsa/*
