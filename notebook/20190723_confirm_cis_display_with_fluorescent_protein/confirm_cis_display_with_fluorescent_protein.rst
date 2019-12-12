********************************************
Confirm CIS-display with fluorescent protein
********************************************

I decided that since my goal is to find out whether or not CIS-display (and 
later: cDNA-display) is working, it would make the most sense to use 
fluorescent proteins fused to repA.  That way I can see the protein directly, 
removing as much ambiguity as possible from the results.

Considerations
==============

Controls
--------
- Reverse-purified fluorescent protein (not fused to repA).  Should be a sharp 
  band.

- Fluorescent protein expressed with unlabeled template DNA.  This would still 
  be a smear, but it make it easier to tell if crosstalk is really a problem.

- FAM-labeled template only.

- FAM-labeled DNA mixed with PURExpress just before loading.  See if PURExpress 
  alone is enough to shift the DNA.

- Stop-codon after mWasabi:

   - See pKBK044
   - Use TAA stop codon: consensus stop codon in E. coli, works with either RF1 
     or RF2.

Fluorescent proteins
--------------------
- Maturation time

   - The gel runs for 2h, so there is some time for the proteins to mature.

- Monomeric

   - My assay is a native gel, so fluorescent proteins that dimerize or 
     tetramerize would affect the results.

      - This could give a bigger gel shift, which would be a good thing.

      - It could also give more smeary bands, which would be a bad thing.

- Brightness

   - My signal is already pretty faint, so I need a reasonably bright protein 
     with a excitation peak that's reasonably well-matched with my lasers.

- Color:
   
   - It's easier to find green proteins (i.e. excited by blue light, 473 nm 
     laser) with good properties.  GFP/Cy5 would be a good pair.

   - My cDNA display linkers have fluorescein, which I can't change.  So if I 
     want something that will work with cDNA-display, I need a protein that's 
     excited by:
     
      - 532 nm (green) light, e.g. a yellow/orange-fluorescent protein.
      - 635 nm (near IR) light

Candidates:

.. datatable:: yfp_candidates.xlsx

Near-infrared candidates:

.. datatable:: nirfp_candidates.xlsx

   Note that these proteins all require the biliverdin cofactor.

Typhoon filters
---------------
I had a hard time tracking down information on the wavelengths permitted by the 
filters available for the Typhoon.  For some filters, I was able to get 
information by searching "<name of filter> filter for fla 9000/9500" on the `GE 
Life Sciences webpage <https://www.gelifesciences.com/en/fi>`_.  For others, I 
tried to interpret the designations in this :download:`Typhoon brochure 
<typhoon_brochure.pdf>`.  The information I compiled is listed below:

Lasers:

.. datatable:: typhoon_lasers.xlsx

Filters:

.. datatable:: typhoon_filters.xlsx

   LP: Longpass filter; DF: Band-pass filter; SP: Shortpass filter.  For 
   longpass and shortpass filters, the number is the cutoff wavelength.  For 
   band-pass filters (I don't know why the abbreviation is DF, but I have found 
   other people using this convention) the numbers are the center and width 
   (FWHM) of the band, respectively.  I don't know what the "R" nomenclature 
   means, but in the case of the BPFR filters, the specific wavelengths 
   transmitted by the filters was available of the GE website.

Azure Sapphire filters
----------------------
The Sapphire has slightly different lasers than the typhoon, but information 
about them was easily available in the control software.  Each laser has only 
one possible filter (this may actually be the case with the Typhoon, too, I'm 
not sure):

.. datatable:: sapphire_lasers_filters.xlsx

   I compiled the list of fluorophores for each laser by looking at page 12 of 
   the Sapphire brochure: :download:`sapphire_brochure.pdf`.  I also compared 
   the imager's lasers and filters with the fluorophore's excitation and 
   emission spectra.  SYPRO Ruby doesn't seem like it'd be excited by a 658 nm 
   laser at all, but if Azure says it'll work, I'll trust them.


Results
=======

Reverse translation
-------------------
The ``reverse_translate.sh`` script automatically includes the proper 5' and 3' 
sequences for Golden Gate assembly into pKBK022::

   ./reverse_translate.sh

Amplify CDS --- 2019/07/24
--------------------------
.. protocol:: 20190724_pcr.txt 20190724_dilute_amplicons.txt

   I accidentally added 2x too much forward primer to the reaction.

   ***

   - Forward primer: ``pUC_seq_amp``
   - Reverse primer: ``ORI_FAM_TM62_REV``, ``REP_FAM_TM62_REV``
   - Template DNA: 23, 24, 25, 26

.. figure:: 20190724_amplify_yfp_repa.svg

- I clearly contaminated my master mix with some template.  I expect that the 
  contamination was one of the templates I used in this reaction (e.g. 23--26).  
  Although these primers would also amplify all of my other repA plasmids, all 
  of the Zif268-repA plasmids would give a noticeably smaller product.  This 
  means that the problem probably isn't that I contaminated my stocks a long 
  time ago, which is good.  I should be more careful next time, though.

- The products should be ~2100 bp (+ORI) and ~1700 bp (-ORI).  I don't know why 
  the products look a little bigger than this.

- The ladder is a little hard to interpret because I ran this gel on top of an 
  old gel, but I think I have the assignments right.

EMSA --- 2019/07/24
-------------------
.. protocol:: 20190724_purexpress.txt

   ***

   Native PAGE:

   - 3-12% gel

   - 2 gels: one with and one without Coomassie in the cathode buffer.

   - Samples:

      - 10 μL IVTT (or ladder)
      - 3.33 μL 4x loading buffer
      - Mix
      - Load 6 μL in each gel

   - Run 150V, 90 min.

.. figure:: 20190724_yfp_repa.svg

- All of the fluorescent proteins emit in the FAM channel (this is a little 
  easier to see by looking at the FAM and HEX channels separately, not shown).  
  This crosstalk is expected for mPapaya and mKO2, but surprising for mApple 
  and mScarlet-I.  I should set up IVTT reactions for mSCarlet-I with and 
  without FAM-labeled template, so I can really see how much crosstalk there 
  is.

- The difficulty of drawing any conclusions with even this small amount of 
  potential crosstalk has convinced me to try this experiment with GFP and Cy5.  
  These two fluorophores are as well-separated as can be, so I should have no 
  concern about crosstalk.  This experiment won't translate directly to my cDNA 
  primers, but I can cross that bridge when I come to it.  Cy5 also isn't quite 
  as bright as FAM, so I might have some trouble seeing the DNA (which is 
  already on the dim side).

- I don't know why there are so many low-MW bands in the HEX channel.  I would 
  imagine that anything in the HEX channel needs to be full-length fluorescent 
  protein.  Maybe these are proteins that didn't finish getting transcribed?  
  That would explain the smears and the low molecular weight.  Maybe I should 
  try a making a C-terminally tagged construct to see if that cleans up the 
  gel.  I'm not sure if a C-terminal tag would interfere with repA, but that's 
  something I could try down the road.

- The differences between the +/- Coomassie gels are interesting:
  
   - Most obviously, the +Coomassie gel is much darker in the HEX channel.  The 
     only explanation I can think of is that the Coomassie is blocking the 
     light that would either be absorbed or emitted by HEX.  This is also 
     consistent with the bottom of the +Coomassie gel being brighter than the 
     top: there's less Coomassie at the bottom of the gel, because most of it 
     gets bound by the ribosomes and stuck at the top.

   - The proteins don't seem to run any differently.  Assuming that the low MW 
     smears in the HEX channel are partially transcribed fluorescent proteins, 
     I wouldn't necessarily expect them to migrate at all in the absense of 
     Coomassie:
     
      - The NativePAGE Bis-Tris buffers have a pH of 6.8.

      - mPapaya and mApple have pI's of about 6.8.  Naively, you'd expect these 
        proteins to not migrate in this buffer, as they should be roughly 
        neutral.

      - mKO2 and mSCarlet-I have pI's of about 5.5.  Naively, you'd expect 
        these proteins to migrate into the gel, as they should be negatively 
        charged in this buffer.

     I can't explain why the two gels aren't more different, especially for 
     mPapaya and mApple.

- The unbound DNA bands are visible in the +Coomassie gel, but not the 
  -Coomassie gel.  The only way I can think to explain this is that the 
  Coomassie is disrupting with the binding between repA and the FAM-labeled 
  DNA.  At the very least, it seems hard to avoid the conclusion that Coomassie 
  may be having some effect on this system, so I should probably stop using it.

- I suspect that the yellow bands about 20% of the way down the high-contrast, 
  -Coomassie gel are the YFP-repA-DNA complex:
  
   - They are yellow, which is what I'd expect the complex to be.  

   - The FAM-labeled DNA has to be in those reactions somewhere.  In panel (c), 
     it's pretty easy to see.  Panel (d) has the same contrast settings---it's 
     part of the same image---so I'd expect it to have about the same amount of 
     green (or more, if Coomassie is interfering with the FAM channel).  Those 
     yellow bands are the only place it could be.

   - Of course, those yellow bands are also pretty bright in the red channel, 
     which makes it hard to discount the possibility of cross-talk.  I'd need 
     the no-FAM control discussed above to eliminate that possibility.
  
   - The corresponding bands in the +Coomassie gel are quite green.  This 
     suggests that these bands really do have DNA, and are not purely 
     crosstalk, especially if my hypothesis about Coomassie blocking the HEX 
     channel is correct.

   - If these bands are the repA complex, it's noteworthy that they run with 
     the ribosomes.  This could be a coincidence: the gel doesn't seem to 
     resolve high-MW species that well, and is clearly overloaded.  But it 
     could also indicate that the complex is not released from the ribosome, or 
     is still being transcribed.  This might be another reason to try a 
     C-terminal tag.

In conclusion, I think I'm seeing binding, but I need to do some more controls 
to be sure.  It'd be more clear with purified protein and a better separated 
fluorophores (e.g. mWasabi and Cy5), but mScarlet-I might be good enough.

EMSA --- 2019/07/26
-------------------
.. protocol:: 20190725_pcr_cloning.txt 20190726_dilute_amplicons.txt 20190726_purexpress.txt

   Printed the wrong protocol by mistake.  Only did the PCR part, not the 
   ligation and transformation.

   ***

   ***

   ***

   Steptactin purification:  

   - Use 1 μL beads per reaction (2 μL total).

   - Follow manufacturer's protocol, except for:

      - Make buffers without EDTA.  Wouldn't have been a problem with 
        fluorescent proteins, but I am planning to do this with Zif268.

      - Elute in 10 μL buffer rather than 25 μL.

   DNase reactions:

   - 5 μL IVTT
   - 1 μL DNase I
   - 0.67 μL 10x DNase I buffer
   - Incubate at 37°C for 10 min.

.. figure:: 20190726_mscarlet_repa_strep.svg

- The high-MW yellow bands are present with or without FAM-labeled DNA, so they 
  definitely represent crosstalk, not binding.

   - This means I definitely need to do my GFP/Cy5 experiment.

- However, I still think this gel shows evidence for binding:
  
   - The template control band is clearly depleted.

   - The low-MW FAM-labeled bands, which are presumably truncated products from 
     the PCR reaction, don't seem to be much affected by the IVTT reaction.  
     This is consistent with these products not containing complete repA/ORI 
     motifs.  In contrast, the full-length template is clearly affected by the 
     IVTT reaction.

   - In is still possible that the full-length template is just being bound by 
     the ribosome or something, and not actually bound by repA.

.. datatable:: 20190726_mscarlet_repa_fam_bands.xlsx

   Gel densiometry results for the FAM channel of the above gel.  Bands #1-#4 
   are just from top to bottom, for those lanes that had multiple bands that 
   seemed to be comparable.

- Note that the +FAM channels do seem to have more FAM fluorescence in every 
  band.  This might imply that the FAM-labeled DNA is distributed through all 
  the peaks.  I wouldn't put too much stock in this, though:
  
   - I didn't quantify it, but the HEX channel also seems brighter for the +FAM 
     lanes, suggesting that the effect is mostly crosstalk.  I don't know why 
     this would be the case, though.

   - The FAM-only control is only 800 px, but the differences between the +/- 
     FAM samples are significantly more than that.  So the differences cannot 
     entirely be attributed to the presence of FAM-labeled DNA.  Most likely 
     the signal is instead crosstalk.

EMSA --- 2019/08/06
-------------------
I had the thought that if the fluorophores are bright enough, I could get 
better data by diluting the reaction so that the lanes wouldn't be so 
overloaded.  Towards this end, I tried serially diluting the IVTT reaction into 
PBS + MgOAc before running the gel.  I used MgOAc because that's what NEB 
recommends for the reverse-His purification protocol, with the comment that the 
Mg helps release product from the ribosome.

.. protocol:: 20190806_pcr.txt 20190806_dilute_amplicons.txt 20190806_purexpress.txt

   ***

   ***

   ***

   - 2x serial dilution of IVTT reactions:

      - Prepare PBS + MgOAc:

         - 100 μL 10x PBS
         - 10 μL 1M MgOAc
         - 890 μL water

      - Put 5 μL PBS + MgOAc in each tube.
      - Transfer 5 μL on each step.

   - Native PAGE:

      - 4-16% gel

      - Each sample:

         - 5 μL IVTT reaction
         - 1.67 μL 4x loading buffer

      - Run 150V, 120 min

.. figure:: 20190806_mwasabi_repa.svg

- The mWasabi-repA fusions don't make it into the 4-16% gel.  I feel like I 
  learned this earlier, then just forgot.  I should add the PAGE parameters to 
  my protocol.

- There is no apparent crosstalk between the GFP and Cy5 channels.  mWasabi 
  also seems to work well.  This should be a good set of fluorophores for 
  testing CIS display.

- I can probably get away with at least a 1/8 dilution.  I want to see what 
  these dilutions look like on a 3-12% gel, though.

- The template DNA is clearly shifted by the IVTT reaction.  This is consistent 
  with CIS-display working as expected, but also with the template DNA not 
  releasing from the ribosome.  There are a number of ways I could 
  differentiate between these two possibilities:

   - Purify the mWasabi-repA fusion from the IVTT reaction.  This is my 
     preferred solution, because a purified reaction would give cleaner data 
     overall.  This is also the only method that would unconditionally tell me 
     whether CIS-display is working.  The other methods would only give firm 
     conclusions if the data looks a certain way.

   - Do a mWasabi (not fused to repA) control.  With no repA domain and no ORI 
     sequence, the DNA should definitely not bind the IVTT product.  If there 
     is binding, it would be attributable to the ribosome (although this 
     doesn't mean that CIS-display isn't also working).  If there is no 
     binding, then CIS_display is working.

     This is basically what the -ORI controls are supposed to do, but I've lost 
     faith in those controls.  The -ORI constructs consistently seem to bind 
     DNA, contrary to my intentions, so they're not good controls because I 
     don't really know how they're supposed to behave.

   - Do the reaction without amino acids.  That would prevent any protein from 
     being translated (particularly the repA domain), but would still allow the 
     ribosome to bing the template DNA.  However, the inability of the ribosome 
     to translate and the absence of release factors may cause the template to 
     bind to the ribosome when it wouldn't normally.

- I should try adding inhibitors of nonspecific binding---BSA, ssDNA, dIdC, 
  Tween---to the IVTT reaction.  These inhibitors are standard components of 
  EMSA protocols, so it would make sense to start using them.  If things are 
  stuck to the ribosome, this might also help unstick them.

EMSA --- 2019/08/07
-------------------
.. protocol::
   
   The same as yesterday, except:

   - 3-12% gel, run for 115 min

.. figure:: 20190807_mwasabi_repa.svg

- The Cy5 didn't seem as bright this time.  That may be because I accidentally 
  left the Cy5-labeled gene on my bench overnight (rather than at -20°C), or it 
  may be that the signal is just more spread out now that it's in the gel.  I 
  increased the contrast in the Cy5 channel to make the bands easier to see.

- The template DNA is clearly gel shifted.  But I don't think this gel tells me 
  anything that the previous gel didn't.

EMSA --- 2019/09/16
-------------------
.. protocol:: 20190916_purexpress.txt

   Used the following genes:

   - no template
   - 27 − T7, i.e. gene amplified with primer 41_TM59 to create a shuffled T7 
     promoter sequence.
   - 44, i.e. mWasabi-repA with a stop codon after mWasabi such that the repA 
     domain isn't translated.
   - 27

.. figure:: 20190916_confirm_cis_display.svg

   −IVTT: Just the DNA templates diluted in water, no IVTT reaction.  +IVTT: 
   PURExpress reactions incubated at 37°C for 2h.  No template: IVTT reaction 
   with no template added.  No incubation: mWasabi-repA (27) template added to 
   a no template reaction immediately before loading the gel.  This controls 
   for the PURExpress buffer conditions (e.g. viscosity, non-specific binding) 
   causing a shift.  No promoter: 27 − T7, shuffled T7 promoter so that no 
   transcription occurs.  mWasabi-STOP-repA: mWasabi-repA template with stop 
   codon after mWasabi.  In this way the template is the same length and the 
   CIS/ORI sequence are still present, but the repA domain is not translated.  
   mWasabi-repA: mWasabi fused to repA.
   
- The DNA templates used in this experiment were designed to be the same 
  length, so that any shift could be attributed to binding.  The −IVTT controls 
  confirm that the templates are in fact the same length.
   
- I think that T7 polymerase is not releasing from the template DNA, most 
  likely because it's bound to the CIS sequence.

  The shift in the DNA is identical between the mWasabi-STOP-repA and 
  mWasabi-repA reactions, so it clearly cannot be attributed to repA binding.  
  The shift is also not present when the T7 promoter has been shuffled, which 
  suggests that the shift is due to the template remaining bound to the T7 
  polymerase.  A reasonable explanation for this is that the polymerase is 
  bound to the CIS sequence, since the role of the CIS sequence is to stall the 
  polymerase [Odegrip2004]_.  I could test this more directly by making 
  constructs with shuffled CIS (and ORI) sequences.

- Reading the methods section from [Odegrip2004]_ again, I also realize that 
  there may be somethings I could do to help release the template DNA from the 
  polymerase:

   - Use the tac promoter (hybrid of trp and lac promoters) with the E. coli 
     RNA polymerase holoenzyme (NEB M0551S).  This is the promoter/RNAP pair 
     that [Odegrip2004]_ used (although they used E. coli lysate, not purified 
     enzyme).  This would test the possibility that T7---which does not work 
     with CIS in nature, unlike the E. coli RNAP---just binds the CIS sequence 
     too tightly.

   - Dilute the IVTT reaction in the same blocking buffer described by 
     [Odegrip2004]_:
      
      - 2-4% (w/v) Marvel: Marvel is a British brand of nonfat dried milk that 
        some people think works the best for blotting.  This brand is not 
        available in HCOM, but I think any nonfat dried milk would be fine.  I 
        could also use BSA, which would probably be even better.
      - 0.1 mg/mL herring sperm DNA: I already have salmon sperm DNA, which is 
        probably fine.  In fact, it seems like herring sperm DNA is considered 
        an older reagent, and salmon sperm DNA is its more modern alternative.
      - 2.5 mg/mL heparin: This could be the important one; heparin occupies 
        the DNA binding site in RNAP, preventing RNAP from binding to promoters 
        and initiating transcription [Wikipedia].  It's conceivable that 
        heparin could also help release RNAP from the CIS-sequence.
      - TBS or PBS

- It's worth noting that in the mWasabi-repA lane, the upper band is green 
  while the lower band is yellow.  In contrast, both bands are about the same 
  intensity of red in the mWasabi-STOP-repA lane.  This suggests that repA is 
  moving the DNA into the lower band.  This phenomenon is in the previous 
  experiment as well.  I'm not sure what the two bands represent, though, so I 
  don't know what to conclude from this.

  Maybe the green band represents extra repA-fusion that's been translated, but 
  that doesn't have a CIS site to bind...

- I still can't firmly conclude whether or not CIS-display is working.  It's 
  promising that the mWasabi-repA fusion superimposes on the template DNA.  But 
  if the DNA isn't released from the RNAP, it's also possible that the ribosome 
  is still bound to the RNAP via the transcript, and that the mWasabi-repA 
  fusion is (for some reason) is still bound to the ribosome.  
  
  To resolve this, I need to release the template DNA/fusion protein from the 
  IVTT machinery.  This was the goal of all the purifications I attempted to 
  do, but they all failed.  I can try the protocols from [Odegrip2004]_ listed 
  above, but other than that I don't have any ideas.
  
- It is possible that some of the protein ran backwards off the gel.  I haven't 
  checked for this possibility in any of my native gels yet.

EMSA --- 2019/10/03
-------------------
.. figure:: 20191003_shuffle_cis_ori.svg

Caveats:

- 55 and 57 both seem shorter than the other constructs.  This is definitely 
  inconsistent with the sequencing data I have for both plasmids, which 
  indicates that the shuffled CIS sequence is the same length as the 
  non-shuffled CIS sequence.  The total GC content is also the same, because 
  the sequence is shuffled (not random).  Sequencing also starts well before 
  the CIS sequence (before oriR, even), and the peaks are very high quality, so 
  this isn't that sequencing missed something.

  I wonder if I used a primer that didn't have the T7 terminator?  I don't 
  think I have such a primer, though.  I'm really not sure what happened, here.

- 56 and 58 didn't amplify well.  This is because I had to design a new reverse 
  primer for the shuffled oriR sequence, and by bad luck the 3' end of that 
  shuffled sequence was *very* AT-rich.  So it's not surprising that the PCR 
  didn't work so well.  I do think these lanes are interpretable, though, 
  because the intended product is clearly present.  I'm assuming that the 
  shorter product is inert.

  Really, the best thing to do about this would be to design a new shuffled 
  oriR with a good primer site at the 3' end.  I could also try gel purifying 
  the right band.  But these results might be good enough.

.. update:: 2019/10/23

   I sent the purified 55-58 plasmids for sequencing, just to make sure I 
   didn't pick the wrong colony to miniprep or something.  The sequencing was 
   noisy across the board, possibly due to contaminants, or posibly just due to 
   bad sequencing.  But the shuffled sequences and their surroundings were 
   entirely correct, so I still don't know what to attribute the apparent size 
   difference of 55+57 to.

Ignoring the caveats discussed above, I can try to draw some conclusions:

- Some part of the transcription/translation machinery binds the oriR sequence.

   - The DNA band is only shifted if both the T7 promoter and the oriR 
     sequences are present.
     
   - I didn't think to shuffle the RBS/start codon, so I can't say if the shift 
     is due to the transcription or translation machinery.  This might be a 
     good thing to test.
     
   - The shift is not due to repA, because it occurs even when repA is not 
     translated.
     
   - The shift is not entirely dependent on the oriR sequence, as there is 
     still a faint shifted band visible even with the shuffled oriR.  This may 
     hint that the shift is due to the transcription machinery, as the 
     translation machinery should disassemblt at the stop codon, which occurs 
     well before. (I confirmed that all of the constructs do have stop 
     codons.)
     
- I think that CIS-display is working, but I still can't be totally sure.
  
   - When the repA domain is not translated, GFP clearly runs independently of 
     the DNA.

   - When the repA domain is translated, a portion of the GFP clearly 
     superimposes on the DNA.  The remainder gets stuck in the well.  Since 
     each DNA molecule can only be bound by one repA molecule, I think this 
     remainder might be excess protein that was translated.  I might be able to 
     get rid of it by running the translation reaction for less time, or by 
     using a greater excess of DNA (I'd rather have excess DNA than excess 
     protein).  Excess DNA might also make the signal between the channels be 
     more even, which would be nice.

   - repA binding (if that is what's happening) doesn't appear to shift the DNA 
     any more than whatever was shifting it already.  

   - Shuffling CIS doesn't seem to affect how much GFP superimposes on the DNA.  
     This is consistent with the fact that repA binds oriR, not CIS.  CIS is 
     believed to stall the RNA polymerase, but that may not be necessary in 
     this system.

   - Even after shuffling oriR, some GFP superimposes on the DNA.  It's not 
     totally clear what's going on here.  A shuffled oriR prevents most of the 
     DNA from shifting initially, and only the DNA that shifts anyways 
     superimposes with GFP.  This could be GFP still associated with the 
     translation machinery, or maybe contamination with DNA that actually has 
     unshuffled oriR?

      - Maybe the brighter yellow bands in the adjacent lanes just mean that 
        more DNA is shifted, so more GFP remains loosely associated with the 
        transcription/translation machinery.

      - The problem remains that I can't distinguish if GFP-repA is "loosely 
        associated with the transcription/translation machinery" or 
        legitimately binding oriR, because I can't get rid of the 
        transcription/translation machinery.  Moving ahead with the qPCR 
        experiments is probably the best way to figure this out.
     
Other observations:

- When neither CIS nor oriR are shuffled but repA is not translated (44), a lot 
  of the DNA gets stuck in the well.  This does not happen for the 
  corresponding construct where repA is translated (27).

- I quantified the intensity of the green fluorescence (mWasabi) of the top two 
  bands for three rightmost lanes:

  .. datatable:: 20191003_shuffle_cis_ori_densiometry.xlsx

     Band 1: Topmost band; basically things that didn't leave the well.  
     Band 2: Second band from the top; presumably things still associated 
     with the transcription/translation machinery.

EMSA --- 2019/11/16
-------------------
I wanted to repeat the above experiment with a shuffled-RBS control to 
distinguish what part of the transcription/translation machinery (e.g. the RNA 
polymerase or the ribosome) is responsible for shifting the DNA.  I also wanted 
to get cleaner PCR products, to make the results a little less ambiguous.

.. protocol:: 20191115_pcr.txt 20191115_dilute_amplicons.txt

   ***

   ***

   .. figure:: 20191115_linearize_shuffled.svg

      All of the amplicons appear to be of the same size, and there was no 
      amplification in the negative control (which had neither template nor 
      primers, so this isn't really saying much).

.. figure:: 20191116_shuffle_rbs.svg

- RNAP (and not the ribosome) appears to be responsible for shifting the DNA.  
  The shuffled RBS control shows a clear shift, which indicates that the 
  ribosome is not the cause of the shift.  Therefore it must be the polymerase, 
  as no shift is seen for the shuffled T7 control.

- RNAP appears to get stalled in the oriR sequence.  This is inconsistent with 
  reports that the CIS element acts to pause the RNAP [Odegrip2004]_ 
  [Praszkier1999]_ [Praszkier2000]_.

  More specifically, [Masai1988]_ details how the CIS sequence contains a 
  Rho-dependent terminator.  The specific requirements of these terminators is 
  unclear, but they are generally understood to comprise an 80-100 bp C-rich 
  region (the "rho utilization site", RUT) followed by a transcriptional pause 
  site.  Rho factor (which is a helicase) binds in the C-rich region of the 
  transcribed ssRNA and eventually displaces RNAP.

  I used [DiSalvo2019]_ to search for putative rho-dependent terminators in the 
  CIS region.  Despite the simplistic nature of this algorithm, I found a 
  reasonable hit that corresponds with the major transcriptional stop site 
  identified by [Masai1988]_.  About 2/3 of this RUT actually occurs in repA, 
  which may explain why shuffling CIS does not have a strong effect.

  Without Rho in the reaction, it makes some sense that maybe RNAP gets stuck 
  on the pause site and fails to release.  The ribosome definitely has a stop 
  codon, so it would be more surprising if it failed to release.  

  It may be prudent to try adding Rho to my IVTT reactions, to see if this 
  helps release my protein from the trancsription/translation machinery.  One 
  way to do this would be to use real cell lysate (rather than PURExpress).  
  Another way would be to add a Rho gene to my PURExpress reactions.  Rho 
  doesn't need any cofactors other than ATP, so this should be enough to get 
  Rho activity.

  All this said, it's still not clear to me why shuffling oriR should release 
  RNAP.
  
  Interestingly, according to Wikipedia__, Rho is blocked by the ribosome but 
  is capable of dislodging polymerase.  In the case of repA, this may be a way 
  to make sure that repA has been completely translated (and therefore given a 
  chance to bind oriR) before the polymerase is dislodged.

  __ https://en.wikipedia.org/wiki/Rho_factor


Heparin incubation --- 2019/09/27
---------------------------------
.. protocol:: 20190927_purexpress.txt

   See binder for PBS-MST recipe, incubation time, native PAGE parameters.

.. figure:: 20190927_incubate_pbs_milk_ssdna_heparin.svg

   MHS: diluted in PBS-MSH (M: skim milk, S: salmon sperm DNA, H: heparin) (+) 
   or just PBS (−).  PURExpress: PURExpress IVTT reaction (+) or just DNA 
   diluted in water (−).

- PBS-MHS did affect the DNA and mWasabi-repA bands, but not in the way I 
  expected.

  It's informative to look at the mWasabi-STOP-repA control.  In this control, 
  the repA domain is not expressed, so the DNA should run to the same place it 
  does in the absence of PURExpress.  However, the DNA is retarded in the 
  PURExpress reaction, consistent with RNAP failing to release from the DNA.  
  Incubating with PBS-MSH retards the DNA even more, so much that it doesn't 
  enter the gel.  This is not consistent with the DNA being released from RNAP, 
  because then it should run the same as it does on its own (I already know 
  that no components of the PURExpress reaction intrinsically retard DNA).  
  However, I also can't conclude that something in PBS-MSH is retarding the DNA 
  (e.g. the milk), the PBS-MSH − PURExpress control doesn't retard the DNA at 
  all.  The only explanation I can think of is that the DNA isn't released, and 
  that PBS-MSH is retarding the RNAP/DNA complex somehow.

- I tried cleaning the glass surface of the gel imager with a microfiber cloth, 
  but it seems to have left little fibers everywhere.  I read that microfiber 
  is supposed to be better than kimwipes, but maybe not.

Heparin incubation --- 2019/09/30
---------------------------------
I'm suspicious that the skim milk is interfering with the gel, so I want to 
repeat the incubation experiment for each buffer component individually.

.. protocol:: 20190927_purexpress.txt

   See binder for PBS-MST recipe, incubation time, native PAGE parameters.

.. figure:: 20190930_incubation_buffer_components.svg

- The smearing in the previous experiment was due to the milk.  Interestingly, 
  this time I see the smearing in both the +/− PURExpress reactions, whereas 
  last time I only saw it in the + PURExpress reactions.  Maybe I did something 
  wrong?

- Heparin and ssDNA both seem to free some DNA from the IVTT machinery (or 
  whatever is running at the top of the gel).  Note the faint bands in the + 
  PURExpress reactions at about the same MW as in the − PURExpress reactions.  
  Since both heparin and ssDNA compete for non-specific DNA binding sites, this 
  freeing effect is consistent with the DNA being bound by T7 RNAP.  Note also 
  that the effect is stronger with both heparin and ssDNA.

  However, the DNA freed by heparin and ssDNA either isn't bound to repA, or 
  doesn't remain bound to repA.  In the mWasabi-repA reactions, all the GFP 
  signal remains is the bands at the top of the gel.  The DNA that is released 
  runs at the same MW as free DNA in the absence of PURExpress, there is no 
  indication it is either bound or retarded by repA.

- Heparin also has an interesting effect on the high molecular weight bands.  
  First, it moves them higher up.  This may be a consequence of the fact that 
  heparin is a polymer and crowding agent.  Second, it causes a significant 
  amount of the mWasabi-repA fusion to never leave the well and to not 
  associate with any DNA.  This could suggest that heparin is also disrupting 
  the repA-DNA interaction.

- Note that 0.8 µL of 75 nM template DNA is 4.5 pmol.  In 
  :expt:`20190626_purify_zif268_repa_via_ribosome_pull_down`, I estimate that 
  each of my 10 µL PURExpress reactions has 24 pmol ribosomes, so there should 
  be an excess of ribosomes.  I'm not sure if there's an excess of polymerase.

- MgOAc has an effect that's similar to heparin and ssDNA, but weaker.  
  According to NEB, MgOAc helps dissociate the ribosome.

I was hoping that this experiment would give me a way to separate the 
repA-complex from the IVTT complex, but unfortunately I still don't see a way 
to distinguish these two possibilities.
