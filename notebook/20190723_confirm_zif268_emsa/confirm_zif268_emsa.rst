*******************
Confirm Zif268 EMSA
*******************

Because I've not seen my Zif268 target being bound by Zif268, I want to confirm 
that I'm doing the gel-shift assay correctly by performing it with Zif268 (not 
fused to repA) purified (via reverse His6) from a PURExpress reaction.

Methods
=======

Clone 35 --- 2019/08/20
-----------------------
I'm having some difficulty with the Golden Gate cloning for pKBK035.

.. protocol:: 20190809_pcr_zif.txt 20190809_pcr_gfp.txt 20190809_golden_gate.txt

None of the colonies I picked had either insert.  Based on the sequencing 
results, I think the sticky ends got filled in by Q5, since I didn't purify the 
PCR products before adding them to the Golden Gate reactions.

.. protocol:: 20190814_golden_gate.txt

   Purify the left-over PCR product from 8/09 and use that for the Golden Gate 
   reaction below.

The colonies I sequenced all had the Zif268 insert but not the mWasabi insert.  
Upon closer inspection, the mWasabi insert seemed to have been replaced by a 
very small insert comprising just the two primers (minus a few nucleotides) 
used to amplify mWasabi.  My only thought is that a primer dimer was amplified, 
then it out-competed the actual mWasabi insert.

.. protocol:: 20190820_pcr_zif.txt 20190820_pcr_gfp.txt 20190820_golden_gate.txt

   ***

   ***

   - PCR cleanup.

   - Run 1% E-gel.

   ***

.. figure:: 20190820_clone_35.svg

   PCR amplified Zif268 and mWasabi inserts.  I increased the constrast to make 
   the low-MW band easier to see.

Both inserts are predominantly correct, but there is a faint low-MW band for 
the mWasabi reaction.  This could be the primer dimer I hypothesized from my 
sequencing results.  Based on this gel, it's probably worth doing a bunch of 
colony PCR reactions to see if I can find one that has the insert.

.. protocol:: 20190821_pcr.txt

   10x colony PCR

.. figure:: 20190821_clone_35_colony_pcr.svg

Most of the colonies seem to have the full length insert.  If so, I guess I'll 
never know what went wrong with the previous reaction.  I picked #1 and #3 to 
grow, miniprep, and sequence.

Reverse-purify --- 2019/08/05
-----------------------------
.. protocol:: 20190805_purexpress.txt

   - SDS-PAGE
      
      - Each sample:

         - 10 μL IVTT aliquot
         - 3.85 μL 4x loading buffer
         - 1.54 μL 10x reducing agent

      - Run 165V, 42 min

.. figure:: 20190806_purify_zif.svg

   SDS-PAGE showing fractions collected during the reverse-His6 purification 
   protocol.  Zif (without any fusions) is 10 kDa.

I don't see any expression of Zif268.  Some thoughts:

- I would expect Zif268 to run about with those dark bands near the bottom, so 
  maybe it's there and just being obscured.  That wouldn't explain why the 
  purification fails, though.  Repeating the experiment with the control 
  plasmid would make this more clear.

- I should make sure my template amplification worked and was clean.

- I should also make sure my template has the right sequence.  The best way to 
  do this would be to make a template plasmid, which I could easily sequence, 
  grow, and use.

- Maybe I should try adding inhibitors of non-specific binding to the reaction 
  buffer, e.g. ssDNA or Tween.  The template DNA doesn't have a Zif binding 
  site, but it may be sticking non-specifically.

Reverse-purify --- 2019/08/19
-----------------------------

.. protocol:: 20190819_purexpress.txt 

.. figure:: 20190820_purify_zif.svg

- This time I used the pKBK034 plasmid as the template.  The advantage of this 
  is that I can be really sure that the sequence is correct.  The disadvantage 
  is that there's probably more RNase in the reaction, so my yields might be 
  worse.

- It's still hard to see, but I think Zif268 is being expressed.  Because 
  Zif268 seems to really overlap with a band in the PURExpress reaction, the 
  negative control is really necessary to see it.

- Next steps:

   - Do a bigger reaction.  It might be that the purification is working fine, 
     and I just can't see it.  This is a 10 μL reaction, diluted to 100 μL for 
     purification.  Maybe I should just bite the bullet and try doing a 50-100 
     μL reaction.

   - Try adding inhibitors of non-specific binding, as described above.

   - Do the gel-shift assay without purification.  I don't think I've tried 
     this with 34 yet.  I tried it with 11 (repA fusion) and it didn't seem to 
     work, but it would still make sense to try it in this context.  I would 
     rather have the mWasabi fusion for this experiment, though.

Reverse-purify --- 2019/08/22
-----------------------------
.. protocol:: 20190822_purexpress.txt
   
   - The below protocol specifies 10 reactions, but really I made a 9x reaction 
     by mixing:

      - 82.8 μL master mix
      - 7.2 μL 75 nM pKBK034
     
     I used the leftover mastermix as a negative control reaction.

   - Not as much of the reaction passed through the spin filter as it usually 
     does (see below).  I would guess that this is because the PURExpress 
     reaction buffer is relatively viscous (although I don't know what's in 
     it).  Perhaps if I do another large scale reaction in the future, I should 
     dilute the reaction 2x no matter what.

      - ~45 μL flow-through
      - ~35 μL retentate

   - When setting up SDS-PAGE, I only used 5 μL of each sample per lane (rather 
     than 10 μL), because these samples are 10x more concentrated.  In 
     retrospect, maybe I should've loaded 10 μL for the flow-through lanes, 
     just because that's where I really wanted to see any small trace of 
     protein.

.. figure:: 20190822_purify_zif_90uL.svg

- The lanes look bad.  Possibly causes:
  
   - Overloaded.  I intentionally loaded a lot more protein than usual, to 
     better pick up faint bands in the flow-through fractions.
     
   - Too much glycerol or something in the loading buffer.  I hardly diluted 
     the reaction into PBSM this time (just 90 μL to 100 μL), and I do think 
     the reaction buffer is fairly viscous.  That could distort the gel a bit.
     
- There seems to be plenty of Zif268 expression.  This agrees with what I've 
  seen previously.

- Most of the Zif268 ends up in the 100K retentate.  Zif268 is definitely small 
  enough to fit through the filter, so this implies that it's binding to 
  something that doesn't go through the filter.  Some possibilities:

   - I confirmed that I ordered the E6800L PURExpress kit, which should have 
     release factors.  So it's not that the protein just isn't releasing from 
     the ribosome.

   - Maybe Zif268 is aggregating.  I haven't seen a visible precipitate, but 
     maybe there's just not enough protein for that.  I don't actually know how 
     stable Zif268 is, and I should double check to make sure I have enough 
     zinc.  Detergent might help with this.  The mWasabi fusion might also 
     help.

   - Maybe Zif268 is binding to a cryptic site in some nucleic acid.  Detergent 
     might also help with this.
  
- The small amount of Zif268 that makes it into the 100K flow-through does not 
  seem to make it into the Ni-NTA flow-through.

   - Zif268 obviously coordinates zinc, and there are reports that it can 
     coordinate nickel as well.  Specifically, both Ni-NTA and Zi-NTA can be 
     used to purify zinc fingers, although different fingers have different 
     affinities for the two metals [Vorackova2011]_  So it may be that Zif268 
     just has intrinsic affinity for the Ni-NTA beads.

   - If Zif268 is just mostly unfolded/aggregated, it may be that whatever got 
     through the filter was just sticking to some other protein component, 
     which itself got pulled down by the Ni-NTA beads.

- It's interesting that the 100K flow-through is usually more-or-less empty.  
  PURExpress has a bunch of His-tagged protein components, so why don't they 
  make it through the filter?  
  
   - Are they all bound in transcription complexes?

EMSA --- 2019/08/28
-------------------
Because reverse-purifying Zif268 has been more challenging than I expected, I 
decided to just do the EMSA experiment with the crude product of an IVTT 
reaction.

.. protocol:: 20190828_purexpress.txt

   See binder for details of binding assay and native PAGE.

.. figure:: 20190828_zif_emsa.svg

- There are faint bands in the Cy5 channel above and below all the target 
  bands.  I would guess that the lower bands are unreacted primers and the 
  upper bands are some kind of mis-amplified PCR product.
  
  If I cared to clean those bands up, I could do PAGE gel extraction.  I don't 
  think an agarose gel would have sufficient resolution.  It might also help to 
  amplify a larger product.  Note that I need to do PCR to attach the Cy5, so I 
  couldn't get cleaner target by cutting it out of a plasmid using restriction 
  enzymes.  I might also try optimizing the PCR.  I already optimized this PCR 
  for yield, but maybe I would get something different if I optimized it for 
  purity.

- I was worried that I might be loading too little material (both DNA and 
  protein) to see, but that wasn't a problem at all.

- The protein ladder was a waste of time.  Without Coomassie, the movement of a 
  protein through a gel is just as much a function of its charge as it is of 
  its mass, so the ladder doesn't really tell me anything.

- The gel-shift effect is very clear, for both Zif268 alone and Zif268 fused to 
  mWasabi.
  
  The shifts cannot be explained by the presence of the PURExpress reaction 
  components---although those components are absent from the "18" 
  lanes---because Zif268 and Zif268-mWasabi evince dramatically different 
  shifts.  Still, it would've been smart to have a "no template" control.

- Despite the large gel shift, the green (GFP) and red (Cy5) bands are not 
  superimposed at all.  This is interesting for several reasons:
  
   - The GFP bands are very heterogeneous.  Why?  Are there lots of incomplete 
     translation products?

   - The DNA bands are also more heterogeneous in the Zif268-mWasabi lanes.  
     (This is easier to see with the green channel turned off.)  This suggests 
     that there are multiple protein species---with slightly different 
     mass/charge properties---binding the DNA.  In contrast, the DNA in the 
     Zif268 lanes appears homogeneous.
     
   - Both Zif268 and Zif268-mWasabi are predicted to be positively charged at 
     pH 6.8 (the pH of the native buffer), while mWasabi is predicted to be 
     negatively charged.  

     .. datatable:: zif_pi.xlsx

         Predicted isoelectric points (pI) from 
         https://web.expasy.org/compute_pi/
     
     This means that both Zif268 and Zif268-mWasabi should migrate upwards off 
     the top of the gel, unless bound to DNA (or something else with negative 
     charge).  What are the green bands, then?  
     
      - They migrate the same in the presence and absence of target, which 
        suggests that they are negatively charged and that they are not 
        necessarily interacting with the target.  Aslo, they migrate faster 
        than the shifted target.

      - They are obviously fluorescent, which suggests that they contain 
        full-length mWasabi.  
        
      - If the bands were mWasabi without Zif268, that would explain the 
        apparent charge and fluorescence, but mWasabi is a C-terminal fusion.  
        So even incomplete translation products should contain the zinc-finger 
        as well.  Maybe there's a cryptic RBS inside the gene?  My DNA 
        templates in these reactions are plasmids, so incorrect PCR products 
        aren't an explanation.

- The addition of BSA and Tween makes no discernible difference for the 
  gel-shift.  On one hand, this isn't surprising.  Both reagents are supposed 
  to help minimize off-target binding, and there are no off-targets in this 
  assay.   On the other hand, some of my other experiments have led me to think 
  that maybe things are getting stuck to the ribosome, and that PBS + Tween 
  might help unstick them.  That doesn't appear to be a problem or a solution 
  here, though.

  Looking at the two "18" bands, though, it's clear that the band with BSA + 
  Tween is much brighter.  Barring the possibility that I just made a pipetting 
  error, maybe BSA or Tween is preventing the DNA from sticking to the plastic 
  tubes.  I remember reading somewhere that DNA tends to adhere to plastic in 
  high-salt solutions, and that a detergent can mitigate the effect.  I can't 
  really explain why, though.
   
EMSA --- 2019/09/18
-------------------
I wanted to repeat the above experiment using both target (TGG) and non-target 
(AAA) DNA, which would better control for different reaction conditions.  I 
also included the Zif268-repA construct.  Previously I had seen that it didn't 
bind it's target, but I want to repeat that experiment with better controls.

.. protocol:: 20190918_purexpress.txt

   Templates:

   - None
   - 34 (plasmid)
   - 11 (gene with RBS)
   - 35 (plasmid)
   - 43 (plasmid)

   My 43 miniprep wasn't good enough to get 75 nM (it was about 50 ng/µL).  
   That isn't a concern for this experiment, though, because that lane is just 
   meant to show where monomeric GFP runs.

   ***

   Binding reaction:

   - 2 µL IVTT
   - 2 µL 6 nM target (same concentration as template)
   - 12 µL PBS-ZBT (see 8/28 for recipe)
   - 8 µL water
   - Incubate at room temperature for 1h.

   Native PAGE:

   - Add 8 µL 4x buffer to each 24 µL reaction.
   - 4-16% gel
   - Load 10 µL/lane
   - Run 150V, 105 min.

.. figure:: 20190918_zif_emsa.svg

   Target: "−" indicates that the middle triplet of the target DNA had the 
   sequence AAA, which is not bound by Zif268.  "+" indicates that the same 
   triplet had the sequence TGG, which is the canonical Zif268 binding site.

- The PURExpress buffers and components themselves do not retard the target 
  DNA.  My previous experiments pointed in the same direction, but this result 
  is much more conclusive.

- Zif268 shifts the target DNA and not the non-target DNA.  Again, the previous 
  experiment demonstrated that Zif268 bound its target in these conditions, but 
  the non-target control makes the result more definitive.

- The Zif268-repA fusion also shifts only the target DNA.  In this case, the 
  target band also becomes noticeably more diffuse than the non-target band, 
  which may reflect that the fusion itself doesn't run as a tight band.

- The protein and target DNA bands are not superimposed.  This indicates that 
  the target DNA is in equilibrium between being bound and unbound, at least 
  until the protein and DNA are separated electrophoretically.

  This can be seen from the Zif268-mWasabi lanes.  The fusion shifts only the 
  target DNA, but the shifted band does not superimpose with either of the two 
  green bands that could represent the fusion.  In fact, the fusion is the 
  green band stuck in the well, because the lower green band runs about the 
  same as the mWasabi monomer, which does not shift the target DNA on its own.  

  I think the most important consequence of this is that I'll need to 
  separately show that Zif268-repA binds its target, and that Zif268-repA binds 
  its encoding DNA.  I've already shown the former here, and I'm pretty close 
  to getting at the latter in my other experiments, so I think I'll be able to 
  move on to my binding assay pretty soon.


Results
=======
- The EMSA experiment seems to be working correctly.

- Zif268 and Zif268-repA are functional.

- Target DNA cannot be used to report on the location of Zif268 in a native 
  gel.

