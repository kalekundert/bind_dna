***********************************
Confirm CIS-display via labeled DNA
***********************************

My goal with this experiment is to observe that CIS-display if successfully 
linking the DNA encoding a protein with that protein.  The approach I am taking 
is to express the Zif268-repA fusion and to run the product on a gel.  I plan 
to visualize the DNA using a fluorescent label (rather than simply staining 
with GelRed).  I plan to visualize the protein by adding target DNA with an 
distinguishable fluorescent label.  Assuming that Zif268 binds and runs with 
its target DNA, both fluorescent bands should be superimposed.

Considerations
==============

Dyes
----
Based on the advice on `this page`__, I decided to used FAM_ and HEX_.  HEX_ 
seems better than JOE_ because its significantly cheaper and doesn't `require 
HPLC purification <https://www.idtdna.com/site/catalog/Modifications/GetAllMods#5>`_.

__ https://www.idtdna.com/pages/education/decoded/article/recommended-dye-combinations-for-multiplex-qpcr

.. _FAM: https://www.idtdna.com/site/Catalog/Modifications/Product/1108
.. _HEX: https://www.idtdna.com/site/Catalog/Modifications/Product/1109
.. _JOE: https://www.idtdna.com/site/Catalog/Modifications/Product/1090

I was considering using Cy3 and Cy5, because those seem to have more separated 
spectra, but they're more expensive and FAM/HEX seem like a pretty common pair.

Ladders
-------
The Invitrogen NativeMark Unstained Protein Standard (for native PAGE) seems to 
include 2 fluorescent proteins, although I can't find a list of exactly what's 
in the ladder.

An alternative ladder I could use is `"Protein Molecular Weight Standards" 
(Serva 39064)`__.  This ladder is also specifically for native PAGE.  It 
contains the following proteins, none of which I believe to be fluorescent:

- Ferritin horse
- Catalase bovine
- Aldolase rabbit
- Albumin bovine (BSA)
- Albumin egg
- Chymotrypsinogen A
- Myoglobin equine
- Cytochrome C

__ https://www.serva.de/enDE/ProductDetails/824_39064_Protein_Molecular_Weight_Standards_html

Troubleshooting
---------------
- Do cDNA display

- Order GFP fusion to repA: confirm CIS display independent of Zif268.  I think 
  this would let me move on to the qPCR experiments with confidence that things 
  are working as they should.

- Make longer target sequence: 

   - Probably the control just ran off the gel, but it'd definitely be better 
     to see it.

   - I'd trust PCR more to give me high-quality dsDNA.

- Try more DNA in reactions.

   - Need a way to see the protein, though, otherwise I won't know if this 
     helps or not.

- Gel gets too hot: proteins denature, run smeary.  Try cold buffers, cold 
  room.  Read more about gel-shift.

- EDTA in gel buffer.  Can't find composition of NativePAGE buffers, but EDTA 
  is a common ingredient.

  .. update:: 2019/07/22

     I asked invitrogen technical support if there was EDTA in the NativePage 
     buffer, and they confirmed that there is *not*:

        Dear Kale,

        Thank you for contacting Thermo Scientific Protein Biology Technical 
        Support.  This is our case 2-5970461952.

        Our NativePAGE running buffer does NOT contains EDTA (or any other 
        chelating agents). I hope this answers your question.

        Please let me know if you have any further questions or comments.

        Sincerely,
        Funmilayo Suleman, Ph.D.
        Technical Applications Scientist III
        Thermo Scientific Protein Biology

     This means I can continue using the NativePAGE gels without worrying about 
     them denaturing my proteins.

     Also, the sepcific composition of the running buffer is described in the 
     NativePAGE manual.  I feel like an idiot for missing that.  It's just 
     BisTris and Tricine.


Results
=======

Amplify target --- 2019/07/19
-----------------------------
I decided to use PCR to amplify and label the target DNA, rather than annealing 
oligos.  This gives me more confidence that I have properly double-stranded 
target DNA.

.. protocol:: 20190719_pcr.txt

   - PCR cleanup, elute in 10 μL EB.
   - Yield: 51 ng/μL

Amplify target --- 2019/07/20
-----------------------------
I was a little surprised by the low yield, so I did an annealing temperature 
titration.  Note the fairly high annealing temperature (which was recommended 
by the NEB Tm Calc) in the previous reaction.  Another possible explanation for 
the low yield is that the silica column didn't retain the product very well.  
The column is advertised as retaining 70 bp to 4 kb, and this product is 68 bp.  
Probably that should be fine, but I'm on the edge.

.. protocol:: 20190720_pcr.txt

   .. note::

      The Eppendorf Mastercycler X50s can do temperature gradients, and it shows 
      clearly which wells are which temperature.  The 96-well Biorad thermocycler 
      can also do gradients, but I wasn't sure which wells were which.  The 
      Invitrogen thermocyclers can't do gradients.

   - Gel densiometry:
      
      - Crop the E-Gel frame.
      - Subtract background: 50 px rolling ball radius.

.. figure:: 20190720_hex_target_optimize_ta.svg

I quantified the intensity of the gel bands to determine the optimal Ta:

.. datatable:: ta_optimization.xlsx

Ta=56°C seems to be optimal for this reaction, and should give at least 2x more 
yield than Ta=64°C.

Control gels --- 2019/07/23
---------------------------
The target didn't run cleanly in either of the two EMSAs I've done, so I wanted 
to run some gels with just a ladder to work out how to run a gel so I can see 
all the bands I want to see.

.. protocol::

   - Lanes:
      
      - 7.5μL NEB 1kb+ DNA ladder, 2.5 μL 4x NativePAGE loading buffer

      - 10 μL NEB Quick-load Purple 1kb DNA ladder (loaded 3 min after the 
        previous lane)

   - Buffers:
      
      - 1x NativePAGE running buffer in both tanks---no Coomassie.

   - 4-16% gel

   - Run gel at 150V until loading dye reached the bottom of the gel (80 min).

   - Stain in 1x PAGE GelRed for 15 min.

.. figure:: 20190723_ladder_only_4_16.svg

- I don't need Coomassie to electrophorese DNA.

- The 4-16% gel doesn't resolve the bands >3 kb very well.  I confirmed that 
  the 1kb+ ladder actually has all the bands I think it should by looking at an 
  old E-Gel.  This actually makes me think I assigned the ladder bands wrong in 
  the 7/22 EMSA.

- I don't know why the 0.5 band for the 1 kb ladder is doubled.  Maybe because 
  it didn't have time to sink to the bottom of the well, since I started the 
  gel immediately.  (Also, the quick-load formulation didn't seem to sink well, 
  as if the formulation was too dense for the small PAGE wells.

.. protocol::

   As above, except:

   - 3-12% gel

   - Only load 5 μL of each ladder + loading dye (i.e. half of above volume)

   - Load target (18) and template (11) DNA:

      - 0.8 μL DNA
      - 1.25 μL 4x NativePAGE loading buffer
      - 2.95 μL water

   - Interestingly, the dye still took about 80 min to migrate to the bottom of 
     the gel.

.. figure:: 20190723_ladder_only_3_12.svg

- The 3-12% gel also doesn't resolve bands >3kb very well.  I'm not really sure 
  why this would be.  My first instinct is that DNA's rod-like shape lets it 
  simply snake through the gel once it gets in the right orientation, and then 
  it's size doesn't matter.  But obviously that's not the case for agarose 
  gels.  It makes me wonder if it wouldn't be better to try using agarose gels.

- The 4-16% gel does give slightly better resolution (e.g. larger distance 
  between 0.1 and 3.0 kb bands), plus I could probably run the 4-12% gel for 90 
  min without the target running off the bottom.  (My guess is that the target 
  would be right at the bottom in 100 min.)  Also, the 4-16% gel is easier to 
  work with.

- This was a good amount of ladder: 3.75 μL ladder + 1.25 μL 4x buffer.

- The 11 CDS runs a little slower than it should.  Perhaps this is due to the 
  FAM modification.

- I really don't know why the 1 kb ladder has two bands at the bottom.  There's 
  only supposed to be 1 0.5 kb band.

EMSA --- 2019/07/10
-------------------
.. figure:: 20190710_zif268_repa.svg

   Native PAGE of reactions with fluorescently labeled template and target DNA.  
   Left panel: Coomassie stain.  Right panel: Composite image showing FAM (red) 
   and HEX (green) fluorescence.  The 4 leftmost (non-ladder) lanes are the 
   PURExpress reactions, and the 6 rightmost lanes are controls containing pure 
   DNA (in the same quantity as added to the PURExpress reaction).  Note that 
   the 11 and 11-ORI templates are 981 kDa and 759 kDa respectively.  Both run 
   faster than the corresponding protein standards, presumably because DNA is 
   globular and perhaps more charged.

- I can see the template DNA (FAM channel), although the bands are faint.

- I cannot see the target DNA (HEX channel), either in the controls or in the 
  IVTT reactions.  The controls are easy to explain: it probably ran off the 
  bottom.  The IVTT results are more complex.  I expected to visualize Zif268 
  with the HEX-coupled target, but that didn't work.  Possible explanations:

   - The target didn't anneal properly.

   - Zif268 didn't bind the target strong enough to keep it from migrating away 
     in the gel.

   - Zif268 isn't folded properly and can't bind its target.

- The IVTT template DNA bands are smeared.  The 11 band is only slightly 
  shifted in the expected direction, while the 11-ORI band is (which I didn't 
  expect to shift at all) shifted significantly relative to its control.

  The shift in the template DNA band in the 11 reaction could indicate that its 
  bound by repA, but it's hard to be confident about that.  Maybe the shift was 
  caused by the ribosomes, somehow.  For example, the ribosomes (or something 
  big) seemed to soak up a lot of the Coomassie, which could cause things to 
  migrate slower.  I think I'll get better results if I can clean up these 
  reactions somehow!
  
- 11-ORI isn't a great control.  I don't really know how it's supposed to 
  behave, and it's consistently not behaved how I would've expected it to.  So 
  I'm not going to put too much stock it it.

- There's something super fluorescent in the ladder!  I found another ladder I 
  could try (see above), but in the future I'll probably leave the ladder out 
  (or maybe run a DNA ladder).

EMSA --- 2019/07/22
-------------------

.. protocol:: 20190722_purexpress.txt

   - Add 0.8 μL 750 nM target DNA to each reaction (10x excess over template 
     DNA)

   - Prepare samples for native PAGE:

      - 10 μL sample
      - 1 μL G-250 additive
      - 3.67 μL 4x loading buffer

   - Native PAGE

      - Use a 3-12% gel.

      - Add 50 μM ZnOAc to the anode and cathode buffers.

      - I tried to load all 14.67 μL of each sample, but it didn't really fit 
        in the wells.  I put a blank lane between each sample to acocunt for 
        this.

      - Run the gel at 150V for 70 min.

   - Image FAM and HEX on the Azure Sapphire

   - Stain in 35 mL 1x PAGE GelRed for 15 min

   - Image GelRed of the Biorad EZdoc.

.. figure:: 20190722_zif_repa.svg

   Ladder: NEB 1kb+.  Note that the high-MW bands are not well resolved, the 
   0.6 kb band is brighter than I'd expect, and the smallest band is largely 
   obscured by the Coomassie.  However, I'm confident in these assignments 
   based on comparisons with the control gels shown earlier.  Left: GelRed.  
   Right: FAM (green), HEX (red).

.. update:: 2019/09/24

   It's interesting to compare this gel with some of my later experiments:

   In this experiment, the target does not appear to be shifted at all.  But in 
   :expt:`35` (which just looks at the gel shift of the target DNA), the same 
   construct (11) seems to shift the target significantly.

   Likewise, in this experiment, the template DNA is shifted slightly and runs 
   as a single band.  But in :expt:`36` (which just looks at the gel shift of 
   the template DNA), the same template is shifted much more dramatically and 
   runs as two bands.

   I think these differences might be due to the use of Coomassie.  Coomassie 
   acts by nonspecifically binding proteins and providing them with a more 
   uniform negative charge, so that they migrate according to their size rather 
   than their innate charge.  However, Coomassie may also disrupt DNA binding.  
   This would explain why the target DNA is not shifted at all in this 
   experiment.  It might also explain why the template DNA is less shifted, 
   although that may also be due to repA having less positive charge in the 
   presence of Coomassie.

   It might be worth doing a side-by-side experiment to know for sure if 
   Coomassie is a problem, but until then, I think I should continue doing all 
   of my EMSA experiments without Coomassie.

- The FAM-labeled template DNA is retarded in both IVTT reactions.  I can think 
  of two ways to explain this:
  
   - repA is actually binding the DNA in both reactions, despite the lack of 
     the CIS-ORI sequence in the 11-ORI reaction.  I don't want to rule this 
     possibility out, because the 11-ORI control hasn't really behaved as I'd 
     expect it to: e.g. not passing through the 100K spin filter.

   - All of the junk in the IVTT reactions is just clogging up the gel and 
     slowing everything down.

  I could test this by adding a lane of IVTT mix plus DNA with no incubation 
  time.  Or maybe IVTT ribosomes but not amino acids or something.

- There is no detectable binding of the HEX-labeled target DNA in either of the 
  IVTT reactions.  I can't say for sure where the Zif268-repA fusion is, but 
  there are no HEX bands other than at the dye front.
  
- I think I should try doing the gel-shift assay with just Zif268 (no repA) 
  purified from an IVTT reaction.  This would more clearly tell me if my gel 
  shift assay is working.

- Although something that looks like the target appears on the Sapphire image, 
  the ladder on the EZdoc image makes it seem like a 70 bp fragment should've 
  run off the bottom of the gel.  The band in the Sapphire image is right at 
  the dye front though; maybe the Coomassie helps the DNA migrate, and the 70 
  bp band is just right at the front?

- I should try running a native gel without Coomassie.  All of the species I 
  can visualize are DNA, and DNA should migrate towards the anode all on its 
  own.  This is just simpler, and it would prevent the excess Coomassie from 
  interfering with the imaging (which I think might be a problem here).  This 
  might even help with the severely overloaded PURExpress lanes, if less of the 
  reaction goes into the gel.  That's probably not so likely, though, because I 
  think most of the overloading is the ribosomes, and RNA is also negatively 
  charged.

- I should also try the 4-16% gel, just the controls, to figure out how long I 
  need to run this gel to see all the relevant DNA bands.


Discussion
==========
I think I tried to test too many things at once with this experiment.  It 
doesn't seem that the target DNA is reliably reporting where the Zif268-repA 
fusion is, and I'm not confident in the 11 - ORI control.  I'm going to try to 
separately test CIS-display and Zif268 binding, in hopes of getting more clear 
answers.

