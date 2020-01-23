************************************
Purify Zif268-repA via affinity tags
************************************

I continue to believe that the PURExpress components are interfering with my 
ability to visualize the success or failure of the CIS-display reactions.  One 
purification method I haven't tried yet is to attach an affinity tag to my 
Zif268-repA fusion.  There are a number of ways I can do this, which I intend 
to explore here.

Note that C-terminal affinity tags would also allow the purification of 
full-length IVTT products from incomplete ones.  This may be especially useful 
for cDNA-display, where the puromycin conjugates could trigger early 
termination.

Considerations
==============

Choice of tag
-------------
- Strep-Tactin:

   - The tag is short, and can be repeated for higher affinity:

      - Strep-tag II: ``WSHPQFEK`` (8 residues)
      - Twin-Strep-tag ``WSHPQFEK-GGGSGGGSGG-SA-WSHPQFEK`` (28 residues)
      - The Twin tag is higher affinity, so I can try that if the shorter 
        tag isn't working well.
      - N- or C- terminal fusion.
      - For the N-terminal fusion, I used the N-terminal sequence from Addgen 
        p29717: ``MAS-WSHPQFEK``

   - Not too expensive, probably about $200 (+ $100 for shipping from Germany) 
     for everything I'd need:

      - $44: Wash buffer (IBA 2-1003-100)
      - $64: Elution buffer (IBA 2-1042-025)
      - $76: Biotin (IBA 2-1016-002)
      - $206: Magnetic beads (IBA 2-4090-002)
      - $112: GFP-Strep-tag control (IBA 2-1006-005)
      - $112: GFP-Twin-Strep-tag control (IBA 2-1007-005)
      - $201: Spin column kit (IBA 2-4151-000)

   - The wash and elution buffers contain EDTA.  Since Zif268 is a 
     metalloprotein, I might want to make these buffers myself without EDTA.  
     The compositions of the buffers are described in protocol, but the only 
     things I'd need would be Tris-HCl (pH=8.0), NaCl, and biotin.

   - The spin columns can elute 50 μL.

   - Assuming that the numbers in the example calculation for determining 
     PURExpress yield are representative, a 10 μL reaction would produce ~100 
     pmol of protein.  The magnetic bead protocol calls for 1 μL beads per 850 
     pmol protein.  So I could use as little as 0.11 μL beads (2.4 μL 5% 
     solution) and elute in as little as 3.0 μL, but I'm sure I'd run into 
     practical problems with volumes so small.  More likely I'd end up using 1 
     μL beads and eluting in 10-25 μL.

   - Nothing from IBA is in the Harvard catalog.  I'm trying to get a quote.

- FLAG:

   - Can be removed by enterokinase.

   - Unlike other epitope tags, FLAG was designed to be convenient for 
     biochemistry (i.e. contains protease cleavage site and is hydrophilic 
     to reduce interference with protein function).

   - c-Myc and HA are other commonly used epitope tags, but FLAG seems to be 
     the most commonly used and best-suited for purification.

   - More expensive.  I'm also a little unclear on exactly what I'd need to 
     buy.  The agarose beads seem to require elution and wash buffers that can 
     only be bought via a larger kit.  The magnetic beads don't seem to require 
     anything but FLAG peptide.  The magnetic beads are also potentially 
     reusable, which is something to think about.

      - $428: Magnetic beads (Millipore M8823-1ML)
      - $481: Agarose beads (Millipore A2220-1ML)
      - $349: 3x-FLAG peptide (Millipore F4799-4MG)
      - $212: FLAG peptide (Millipore F3290-4MG)

- MBP:

   - Not great to have a big tag for IVTT.
   - Clone with TEV site so removal is possible.
   - Not so excited about this possibility, but it would be a nice contrast 
     to all the short tags I'm thinking about.

- His:

   - This would conflict with the components of the PURExpress reaction, but it 
     would be possible to purify my proteins by relying on the fact that they 
     shouldn't pass through the 100K spin filter.

   - Although the 11- ORI control seems weird.


Methods
=======

EMSA --- 2019/07/26
-------------------
See :expt:`20190723_confirm_cis_display_with_fluorescent_protein`.

.. figure:: 20190726_mscarlet_repa_strep.svg

The streptactin purification didn't work at all.  Some things I could try:
  
- N- and C-terminal tags.
- Twin tags
- Columns rather than batch.  Would be annoying to have to get another 
  quote, but that's just how it is.

Streptactin --- 2019/08/21
--------------------------
.. protocol:: 20190821_dilute_amplicons.txt 20190821_purexpress.txt

   ***

   ***

   Streptactin purification: 

   - Based on 7/26 protocol.

   - For each reaction:

      - Wash 2 μL beads (40 μL 5% bead solution).  I washed the beads 
        separately this time, for no particular reason.  Keep washed beads on 
        ice.

      - Dilute the reaction to 50 μL with wash buffer - EDTA (i.e. add 40 μL 
        wash buffer - EDTA).

      - Save 10 μL "crude" aliquot.

      - Add remaining 40 μL to washed beads.  Vortex to resuspend.

      - Incubate on ice for 30 min.  Vortex every 5 min.

      - Wash three times with wash buffer - EDTA:

         - First wash: 50 μL buffer, save 10 μL "wash 1" aliquot.
         - Second wash: 100 μL buffer, add 10 μL to "pooled wash 2,3" aliquot.
         - Third wash: 200 μL buffer, add 10 μL to "pooled wash 2,3" aliquot.

      - Elute with 50 μL eultion buffer - EDTA.

   Electrophoresis:

   - Run SDS-PAGE with aliquots from purification.

   - Run native PAGE with the eluted material.

.. figure:: 20190821_strep_tags.svg

- Interestingly, the Twin-Strep tag fusions run noticeably slower than the 
  single Strep-tag fusions (gel A).  I didn't expect such an easily observable 
  difference, but the molecular weights seems just about right.

- All of the expressed protein seems to be lost in the initial flow-through 
  (gel A).  No product is visible in any of the washes (gel B) or in the 
  elution (gel A).  This indicates than none of the protein binds the beads.

- Neither putting the tag on the N- or C-terminus nor using the "single" or 
  "twin" tag seems to improve the purification.

- The native PAGE gel (gel C) is interesting.  None of the eluate lanes have 
  any signal in the Cy5 channel, meaning that no repA-gene complex was 
  purified.  However, it's notable that the lanes with N-terminal tags (and not 
  the C-terminal tags) contain a number of faint green bands.  I think these 
  bands represent partially translated products: long enough to be fluorescent, 
  but not full-length (otherwise there wouldn't be multiple bands).  
  
  In fact, looking at these bands, I can see *very faint* corresponding bands 
  in the SDS-PAGE eluate lanes (gel A), around 10 kDa.  I don't have a protein 
  ladder in the native PAGE gel (and such a ladder might not correspond to 
  actual size very well anyway), so I can't really say how big those green 
  bands are.  For reference, mWasabi is 27 kDa (the mWasabi chromophore is 
  about 1/4 way into the sequence) and repA is 33 kDa.  So if the green bands 
  are really ~10 kDa, they would not even be full-length mWasabi, but they 
  would at least contain the chromophore.

  If I can purify partially expressed protein, that suggests that there's 
  something about repA that is interfering with the purification.  It might be 
  worth cloning Strep Tag + mWasabi without repA to probe this more directly.  
  It's also worth thinking about what about repA could be interfere with 
  binding, and what I could do about it.  Some things that come to mind:
  
  - Non-specific binding to DNA or the ribosomes that blocks access by 
    streptactin.
    
  - The negative charge of the bound DNA repels streptactin or adheres to the 
    tag.

Streptactin --- 2019/08/30
--------------------------
Address the questions raised on 8/21 by purifying Dual Strep-mWasabi with and 
without the repA C-terminal fusion.

.. protocol:: 20190830_purexpress.txt

   Used 38 and 38-repA for this experiment (N-terminal dual Strep-tagged).

   ***

   - Streptactin purification (see 8/21)

      - Wash 2 μL beads (40 μL bead solution)

      - For each reaction:

         - Dilute to 50 μL with wash buffer

         - Save 10 μL aliquot ("crude")

         - Add reaction (40 μL) to 1 μL washed beads.

         - Incubate on ice for 30 min.  Flick to mix every 5 min.

         - Save 10 μL aliquot ("flow-thru")

         - Wash with 40 μL wash buffer

         - Save 10 μL aliquot ("wash 1")

         - Wash again with 40 μL wash buffer.

         - Save 10 μL aliquot ("wash 2")

         - Add 25 μL elution buffer.

         - Incubate 10 min.

         - Remove beads.

         - Save 10 μL aliquot ("eluate")

   - SDS-PAGE

      - Samples:

         - 10 μL IVTT
         - 3.85 μL 4x buffer
         - 1.54 μL 10 reducing agent
         - 70°C for 10 min.

   - Native PAGE

      - DNA:
         
         - 0.8 μL 75 nM DNA
         - 49.2 μL water

      - Samples:

         - 10 μL IVTT/DNA
         - 3.3 μL 4x sample buffer

.. figure:: 20190830_mwasabi_strep_tag.svg

   (a) Coomassie gels of fractions collected during the Streptactin 
   purification.  (b, c) Native PAGE gels of the +/- repA eluates.  Both panels 
   are the same image, just with different color balance settings applied.

- mWasabi (without the repA fusion) is expressed and purified well.  This 
  indicates that the Streptactin beads are functional.

- mWasabi-repA is expressed well, but is not bound by the beads.  Nothing is 
  purified and all of the protein can be seen in the flow-through fraction.

- The mWasabi-repA lanes look about the same as in the previous experiment.  
  The faint green bands do seem to correspond to incomplete translation 
  fragments that contain most/all of mWasabi and little/none of repA.  This 
  supports the hypotheses that something about the full length fusion protein 
  interferes with Streptactin binding.  I still don't know the exact reason for 
  the interference, but I think I can conclude that I will not be able to 
  purify the repA fusion using Streptactin.
  
  Note that how a protein runs in a native gel depends on its size and charge, 
  so I should be careful to avoid concluding too much from gel shifts that 
  could be caused by multiple factors.

His6-TEV --- 2019/09/10
-----------------------
Purify GFP-repA with Ni-NTA

.. protocol:: 20190910_pcr.txt 20190910_dilute_amplicons.txt 20190911_purexpress.txt

   More details in binder.

   Ni-NTA purification:

   - Dilute rxns to 250 μL with lysis buffer.

   - Aliquot 10 μL (crude)

   - Add 50 μL magnetic Ni-NTA beads.

   - Incubate 1h, 4°C, continuous mixing

   - Remove beads

   - Aliquot 10 μL (flow-thru)

   - Wash with 250 μL wash buffer

   - Aliquot 10 μL (wash 1)

   - Repeat above wash

   - Aliquot 10 μL (wash 2)

   - Add 25 μL elution buffer

   - Incubate 2 min, 4°C

   - Aliquot 10 μL (eluate)

   - Setup TEV reaction

      - 10 μL eluate
      - 1 μL 0.5 mg/mL EZCut TEV protease

   - Incubate 1h, 34°C

   - Dilute to 125 μL with PBST

      - PBST being the base for the Qiagen Ni-NTA purification buffer without 
        any imidazole.

      - This brings the final imidazole concentration to 20 mM, same as in the 
        wash buffer.  Might be better if it were 10 mM (same as the lysis 
        buffer), but if I dilute the reaction too much I won't see anything.

   - Add 50 μL magnetic Ni-NTA beads

   - Incubate 1h, 4°C, continuous mixing

   - Aliquot 10 μL (tev)

.. figure:: 20190911_purify_mwasabi_repa_via_his6.svg

   41-repA: His6-TEV-mWasabi (e.g. just mWasabi); 41: His6-TEV-mWasabi-repA 
   (e.g. N-terminal tag); 42: mWasabi-repA-TEV-His6 (e.g. C-terminal tag).  (a) 
   Given protein MWs (left) include the His6+TEV tag.  (b) DNA: Only 
   Cy5-labeled DNA; eluate: Ni-NTA eluate from purification.
   
- The 41-repA band runs slightly slower than I'd expect based on the ladder.  
  This band is faintly visible in the crude reaction, and clearly visible in 
  the eluate.

- The 41 and 42 bands are visible in the crude reaction and the flow-through 
  fraction, indicating that the repA fusions did not bind the beads.  This is 
  supported by the fact that no high-MW mWasabi band is see in the native gel.

- There is are bands in the 41 and 42 eluate lanes that line up pretty well 
  with the mWasabi-repA fusion, but these bands are also present in the two 
  control reactions, so I do not think it actually represents the fusion.

- Some low-MW species are purified from both reactions with N-terminal tags.  I 
  assume that these species are basically just N-terminal fragments that got 
  created due to some error.

- The purification works well on mWasabi without repA, suggesting that 
  something intrinsic to the repA fusions is making them hard to purify.  This 
  is the same result I got with Streptactin.
  
- Based on the native gel, it looks like I can purify small amounts of what 
  seem to be partially translated products with the N-terminal tag.  This is 
  also consistent with the idea that the repA fusions interfere with 
  purification, since I do not see full length protein.

- The TEV cleavage didn't seem to work, but it might just be that I diluted the 
  reaction too much to see anything.

TBS-HMST --- 2019/10/01
-----------------------
After getting some evidence that heparin, ssDNA, and MgOAc may help release the 
repA-DNA complex from the IVTT machinery (see 
:expt:`20190723_confirm_cis_display_with_fluorescent_protein`), I want to try 
purifying the repA complex after incubation in this buffer.  This experiment is 
a bit of a crapshoot.  If I am able to purify a repA-DNA complex, that would 
provide solid evidence both for CIS-display working and for the incubation 
buffer helping.  If I'm not able to purify the complex, I won't really learn 
anything.

.. protocol:: 20191001_purexpress.txt 20191001_strep_tactin_purification_with_aliquots.txt

   Templates: 27, 38, 39, 40

   See binder for TBS-HMST recipe.  I didn't use BSA ("B") because I didn't 
   want to see it on the gel, and it didn't seem to be important from previous 
   experiments.  I used Tris + NaCl (TBS) rather than PBS to match the buffers 
   used for Strep-tactin purification.

   See binder for purification protocol.  I made diluted reactions to 100 µL 
   (as in my previous incubation experiments) rather than 50 µL (as in my 
   previous purificaiton experiments).  To compensate for the more dilute 
   samples, I took 15 µL aliquots (rather than 10 µL).

   See binder for SDS-PAGE protocol.  Because I diluted my samples more, I 
   loaded 20 µL, which is more than usual.

.. figure:: 20191002_purify_strep_mwasabi_repa_with_pbs_hmst.svg

   Lanes: 27: N-terminal Strep-tag.  38: N-terminal Twin Strep-tag.  39: 
   C-terminal Strep-tag. 40: C-terminal Twin Strep-tag.

- I can't say much about this data because I don't see the bands for the 
  mWasabi-repA fusion anywhere.  The gel is faint because everything is more 
  diluted than usual, and I forgot to take aliquots of the crude reaction, so I 
  don't really know if I even should be able to see the repA fusion.  If the 
  protein were visible, I'd expect to see it around 62 kDa, as I've seen in 
  previous purifications (e.g. 8/21).

  I also didn't include the STOP control, which would've verified that the 
  beads were working.  I could've set up this experiment much better.

- The repA fusion is certainly not visible in the elution lanes.  This isn't 
  really attributable to dilution, because the protein is eluted in 50 µL.

Although this wasn't a great experiment, there's no evidence that the 
purification worked.


S30 extract --- 2020/01/21
--------------------------
As described in :expt:`20190723_confirm_cis_display_with_fluorescent_protein`, 
I think that rho factor may help disassemble the transcription/translation 
machinery.  Since my current hypothesis is that the affinity tag is being 
sequestered by said machinery, expressing the protein in S30 extract (which has 
rho factor) might allow the purification to succeed.

.. todo::

   Purify with S30 extract.

Results
=======
It does not seem possible to purify repA fusions using purification tags.  I 
don't know what the problem is, but my guess is that is has something to due 
with repA preventing the transcription/translation machinery from dissembling 
and just sterically prevent the tags from being bound.

Going forward, I don't really need to purify the repA constructs in order to 
use them.  Instead, I'll need to run an experiment (already planned) with some 
really good controls to show that repA is binding DNA.  From there, I can still 
try to use repA in my binding assay.

It's worth noting here that I also didn't observe any DNA-binding by Zif268 
when fused to repA.  I should double-check that experiment, but it may be that 
repA fusions are just really not functional for some reason.

I also want to try purifying cDNA my constructs, since I anticipate that those 
will be better behaved.
