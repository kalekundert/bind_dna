********************
Confirm cDNA display
********************

I want to see whether or not I can successfully create a linkage between attach 
a translated protein to the cDNA encoding it.  This is the first step in doing 
cDNA-display-based DNA-binding assays.  My goal is to run a gel where both the 
protein and the cDNA can be seen in different channels, and to show that the 
two bands are superimposed.  As a negative control, I can express the protein 
without the linker.  I'll mostly be following the cDNA-display protocol from 
[Naimudden2016]_, with some PURExpress-specific details taken from 
[Barendt2013]_.

Considerations
==============

Visualization
-------------
To see whether or not CIS-display was working, I expressed GFP using DNA 
labeled with Cy5, so I could simultaneously visualize both molecules.  See 
:expt:`20190723_confirm_cis_display_with_fluorescent_protein` for details.  To 
see whether or not cDNA display is working, I'll have to do things a little 
differently, because the linker I ordered is labeled with FITC.  This means 
that I have to visualize the protein in the Cy5 channel.  I thought of a few 
ways to do this.  Of these ideas, I think the fluorescent protein stain seems 
the most promising:

- Use a fluorescent protein stain.
  
   - Because cDNA display creates a covalent attachment, I should be able to 
     use SDS-PAGE to see whether or not the attachment was formed.

   - Because SDS-PAGE is so much higher resolution than native PAGE, I will 
     probably be able to see a distinct band for my product even if I just 
     stain all proteins.  This would remove the need to do any purification 
     steps (although such steps would still give nicer data).

   - Because I wouldn't have to use a fluorescent protein, I could just use my 
     Zn-finger, which would make my result more relevant to the actual binding 
     assay I want to perform.

   - While I could use Coomassie, a fluorescent stain will be more sensitive, 
     and will allow me to create nice multichannel images using the gel 
     imager.
     
   - A :download:`brochure for the Sapphire imager <sapphire_brochure.pdf>` 
     claims that SYPRO-Ruby can be visualized using the Cy5 laser.  I'm a 
     little skeptical because SYPRO-Ruby doesn't appear to be excited by 658 nm 
     light, but I trust the brochure enough to give it a try.
     
   - SYPRO-Ruby might overlap with the GFP channel.  It will be almost 
     maximally excited by the 488 nm laser, but it shouldn't emit in the 
     507-529 nm band.  If there is cross-talk, I might also try SYPRO-Ruby, or 
     maybe the Typhoon will have better laser/filter combinations.

- Use a near-IR fluorescent protein:

   - I don't want to use an RFP, because most RFPs have some overlap with the 
     GFP channel.  Even if the overlap is small, it makes the results very 
     difficult to interpret.

   - Near-IR proteins all seem to require the biliverdin cofactor.  This could 
     cause issues with background fluorescence or getting the protein to fold, 
     although most likely it would be fine.

   - I would need to continue using native gels, so that the fluorescent 
     protein would remain folded.

- Western blot

   - Kinda the best of both worlds, in that I can run an SDS-PAGE gel and then 
     visualize just my protein of interest.

   - I'd probably just use an epitope tag to find my protein, which I could 
     possibly also use for purification (see below).

   - I haven't done a lot of Westerns before, so I'd have to learn how to do 
     that.  Could be a bit of a learning curve, but it'd also be a good 
     technique to know.

- Order a new linker with Cy5 instead of FITC.

   - This would be expensive.

   - But it would make it possible to visualize the mRNA/linker conjugation 
     using SYBR gold.

   - But then I couldn't use SYPRO-Ruby, which could be really nice.

   - If I were to do this, though, I should also think about ordering the 5' 
     end phosphorlyated.  That would let me skip the phosphorylation step, and 
     would only improve efficiency.

Purification tags
-----------------
Since I'm planning to do SDS-PAGE, it shouldn't be absolutely necessary to 
purify my proteins.  But I am worried about getting partially transcribed 
products due to puromycin being present in the reaction.  To account for this, 
I'd like to have a C-terminal tag I can use to purify fully-translated 
products.  

Note that this C-terminal tag will really be an internal tag, because the 
puromycin linker will be attached to the C-terminus.  This limits the selection 
of tags I can use.  I consulted the list of affinity tags compiled by 
[Kimple2013]_ to find the best candidates.

- His-tag: [Naimudden2016]_ used a C-terminal His-tag (flanked by GGGS motifs 
  on either side) and seemed to be able to successfully purify protein, even 
  though [Kimple2013]_ claims that His-tag are not internal tags (although I 
  don't see why they would need to be at either terminus to coordinate Ni).  Of 
  course, all of the proteins in the PURExpress reaction are also His-tagged, 
  but I could always use a 100 kDa spin-filter get rid of most of them 
  (although not any that might run near my protein in a gel).

- Strep-tag: I don't really have any reason to believe that this would work 
  internally, except that the Twin Strep-tag has one Strep-tag that's basically 
  internal, and it still seems to contribute something to binding.  Plus I 
  already have all the beads and buffers, so I might as well try it.

- HA-tag, Myc-tag: These are two common epitope tags that are considered by 
  [Kimple2013]_ to be internal tags.  Antibody purifications are considered to 
  give low yields, but I think they're worth trying if the above tags don't 
  work.

C-terminus
----------
[Naimudden2016]_ calls for the following sequence elements C-terminal of the 
gene of interest:

- 2x GGGS
- His-tag
- 1x GGGS
- Y-tag

I'm curious if these C-terminal elements can't be improved.

- Stop codon with ΔRF PURExpress

   - The idea with using a stop codon is that, if there are no release factors 
     in the reaction, the ribosome will stall on the stop codon, which will the 
     puromycin a chance to react.

   - I might even be able to subsequently add release factors to help release 
     the puromycin-linked protein from the ribosome.  None of the papers I've 
     read have seemed to say much about releasing from the ribosome, so maybe I 
     don't need to worry about it.  But if it turns out to be a problem, this 
     could be a solution.

   - Using a stop codon would require me to use PURExpress, because any lysate 
     would include release factors.

   - A stop codon would make the C-terminus more consistent.  Currently I 
     believe that translation will stop when the ribosome encounters the 
     double-stranded region where Linker-N is bound, but it's also possible 
     that Linker-N itself will be translated (it's mostly Gly).  A stop codon 
     would make it more clear where translation stops.  
     
     Unfortunately, though, I don't really have a way of knowing how consistent 
     the C-terminus is.  So I won't know if this is really a problem, or if a 
     stop codon even does anything to help.  If I really wanted to know this, 
     I'd probably have to do mass spec.

   - A stop codon would also enable the use of 3' barcodes.  See 
     :expt:`20190403_detect_binding_via_cdna_display`.

   - The poly-A arm supporting the puromycin in Linker-N is optimized for use 
     without a stop codon.  Granted, it's also optimized for use with 
     eukaryotic ribosomes, which are slightly bigger than bacterial ones.  
     Still, adding a stop codon would move the puromycin further from the 
     A-site, which may interfere with the ability of puromycin to attack the 
     stalled peptide.

  All told, using a stop codon could be a good thing or a bad thing.  It's not 
  something that [Naimudden2016]_ would've tried, because they were using 
  rabbit reticulocyte lysate.  So I think it's worth trying, but I should get 
  the system working first without it.

- Remove or shorten the GGGS linkers

   - First of all, [Naimudden2016]_ used GGGS linkers, while most everyone else 
     uses GGGGS (4 Gly instead of 3).  It probably doesn't matter, but I'm 
     tempted to use the more canonical (longer) linkers.

   - I'm skeptical that 2 repeats of the Gly/Ser linker are needed before the 
     His-tag.  After all, the His-tag and the second linker can be thought of 
     as one big linker between the protein and the mRNA/cDNA.  Again, though, 
     this probably doesn't matter too much.

Variable transcript length
--------------------------
As described very clearly by [Gholamalipour2018]_, T7 polymerase can add a 
variable number of extra nucleotides to the end of transcripts.  This is caused 
by the transcripts folding on themselves and priming additional transcription.   
This phenomenon occurs more readily when the concentration of RNA is high (as 
in IVT reactions), and in these conditions [Gholamalipour2018]_ show that the 
correct-length transcripts can become a minor product.

It's possible that this extra transcription could interfere with the ligating 
of Linker-N to the 3' end of the transcript.  However, the actual Y-tag 
sequence seems like it would be fairly resilient to this effect, because it it 
mostly G and therefore has poor self-complementarity.  Also, [Naimudden2016]_ 
reported good ligation yields, and they didn't mention this issue at all.  So 
I'm not going to worry about this for now, but I wanted to at least put my 
thoughts down.

Linearization
-------------
[Naimudden2016]_ linearized their DNA template by using PCR to add the 
promoter, the GS-linkers, the His-tag, and the Y-tag to the coding region.  I 
think that makes sense for directed evolution, where you'll only get out the 
coding sequence after each round.

In my case, I designed my plasmids to include all of the aforementioned 
elements.  I think this approach makes sense for me, but there are some 
trade-offs.  The advantage is that I can verify the sequence of all of those 
elements in advance, and I don't need to do PCR with ridiculously long primers.  
The reverse primer [Naimudden2016]_ used must have been >100 nt long, and the 
forward primer wouldn't have been short either.  The disadvantage is that I 
can't exactly control the 3' sequence.  I tried using PCR to amplify the 
cassette out of the plasmid, but the reaction didn't work well because the 
Y-Tag is almost 80% GC.  The GS-linkers are also quite GC-rich, but I could try 
making a ~60 nt primer that binds in the His-tag.  Alternatively, I put an XmnI 
site just after the Y-Tag.  This robustly linearizes the plasmid, but leaves an 
extra nucleotide.  I don't think there are any enzymes that can make a blunt 
cut at the right position, although it might be possible to use a Type-IIS 
enzyme to make an sticky-end cut and to subsequently either cleave off or fill 
in the overhang.  I don't think the extra nucleotide is going to cause any 
problems, but it's something to think about if I need to troubleshoot.

Ligation
--------
The protocol given by [Naimudden2016]_ for ligating the linker to the 
transcribed RNA has some ambiguities:

- The buffer for resuspending the linker is not specified.  EB is probably 
  reasonable, though.

   - The "Certificate of Analysis" from Midlands CRC (which is saved in the 
     binder) specifies that the oligos should be resuspended in a "sterile 
     buffer with pH ranging from 5--9" and stored at -20°C.  So EB will be 
     fine.

- The buffer for the annealing and ligation reactions is not specified.

   - T4 PNK (NEB M0201) requires 1x T4 *DNA ligase* buffer:
     
      - 50 mM Tris-HCl
      - 10 mM MgCl₂
      - 10 mM DTT
      - 1 mM ATP
      - pH 7.5 @ 25°C
        
   - T4 PNK (NEB M0201) can also be used with T4 PNK buffer, provided that ATP 
     is also added:
     
      - 70 mM Tris-HCl
      - 10 mM MgCl₂
      - 5 mM DTT
      - pH 7.6 @ 25°C

   - T4 RNA ligase (Takara 2050A) conveniently requires almost the same buffer 
     as T4 PNK:
     
      - 50 mM Tris-HCl
      - 10 mM MgCl₂
      - 10 mM DTT
      - 1 mM ATP
      - 0.01% BSA (recommended)
      - pH 7.5 @ 25°C

   - Considerations in choosing which buffer to use:

      - Salt is needed for annealing, but is inhibitory for phosphorylation.

         - The annealing protocols I found (discussed below) all call for 
           50--150 mM salt, typically NaCl.  This makes sense to me, because 
           the salt presumably helps the two polar molecules come together.

         - 150 mM NaCl inhibits 50% of T4 PNK activity:
           https://www.neb.com/faqs/2011/11/22/what-factors-can-cause-incomplete-phosphorylation-when-using-t4-polynucleotide-kinase

         - T4 RNA ligase is inhibited by metal chelators, so I should not use any 
           buffers with EDTA.  But NaCl is not mentioned as an inhibitor.

      - ATP

         - At first, I was worried that if I added ATP before the annealing 
           step, it would break down during the extended high-temperature 
           incubation.  But it occurs to me that the NTPs in PCR are stable 
           after similar incubations at similar temperatures, so probably I 
           don't need to worry about this.

   - Strategies:

      - Phosphorylate the DNA linker in advance, then add as much salt as I 
        want to the annealing reaction.

      - Drop dialysis after annealing to remove salt.

      - Do the annealing reaction in a small volume (e.g. 5 µL), then dilute 
        for the ligation reaction (e.g. to 50 µL) to reduce the salt to a 
        non-inhibitory level.

   .. update:: 2019/12/11

      [Naimudden2011]_ clarifies that "mRNA was annealed to the biotinylated 
      puromycin-linker DNA (1:1 ratio) via the Y-tag sequence in 1× ligase 
      buffer (Takara, Kyoto, Japan)..."  So they added ATP before the annealing 
      gradient, and didn't add any salt for the annealing.  It seemed to work 
      fine for them though, so maybe I can just do this moving forward.

- The volume of the annealing and ligation reactions is subtly implied, I 
  think.  It is specified that 50 pmol of mRNA and linker-N are added to the 
  reaction.  It is later specified that the mRNA/linker conjugates are added to 
  the IVTT reaction at a concentration of 100 nM.  50 pmol at 100 nM would 
  correspond to a volume of 500 µL.  That's probably not right, though.  
  Perhaps 100 nM is a final concentration, or not all of the conjugate is added 
  to the IVT reaction.

  [Naimudden2011]_ has a little more detail.  First, the mRNA/linker conjugate 
  is purified using an RNeasy kit after ligation.  Then 3-5 pmol of the 
  conjugate are translated in a 25 µL retic lysate reaction.  The manual for 
  that kit calls for up to 5.75 µL of RNA per 25 µL reaction, so after 
  purification the concentration is about 1 pmol/µL (1 µM).

  Because of the RNeasy step, this doesn't say anything about the volume of the 
  ligation reaction.  Annealing supposed works best with high concentrations of 
  oligos.  Ligation might be better at lower concentrations, to avoid ligations 
  between mRNA molecules.  If that's a problem, though, I could avoid it 
  entirely by phosphorylating the DNA linker and dephosphorylating the mRNA.  
  Or in the future, ordering a phosphorylated linker.  For now, I'll just do 
  something reasonable like 10 µL.

- "mRNAs were annealed to linker-N by heating at 94°C and gradient-cooling to 
  4°C."

   - The duration of the gradient cooling is not specified.

   - `This protocol 
     <https://www.sigmaaldrich.com/technical-documents/protocols/biology/annealing-oligos.html>` 
     from Sigma specifies:
      
      - Thermocycler protocol:

         - 95°C for 2 min
         - 95°C→25°C over 45 min.
         - Hold at 4°C

      - Buffer:

        - 10 mM Tris
        - 50 mM NaCl
        - 1 mM EDTA

      - Oligo concentration: 50 µM as an example, didn't seem like strict 
        requirement.

   - `This protocol 
     <https://www.idtdna.com/pages/education/decoded/article/annealing-oligonucleotides>` 
     from IDT specifies:

      - 94°C for 2 min, the "gradually" cool.  Specific cooling time not given, 
        but mentions that you can just take the reaction out of a heat block 
        and leave on the bench.

      - Buffer: "This provides a buffering environment and the salt is 
        necessary for oligonucleotide hybridization."  Available for purchase 
        from IDT.

         - 100 mM KOAc
         - 30 mM HEPES
         - pH 7.5

      - Oligo concentration: 10--100 µM

   - `This protocol 
     <https://tools.thermofisher.com/content/sfs/brochures/TR0045-Anneal-oligos.pdf>` 
     from Thermo specifies:

      - Thermocycler protocol:

         - 95°C for 5 min
         - 95°C→25°C over 70 min (1°C/min)
         - Hold at 4°C

      - Alternative thermocycler protocol; pause at annealing temperatuer (Ta): 

         - 95°C for 5 min
         - 95°C→Ta at 1°C/min
         - Ta for 30 min
         - Ta→25°C at 1°C/min
         - Hold at 4°C

      - Buffer: "Tris or phosphate buffer containing salt; for example, 10 mM 
        Tris, 1 mM EDTA, 50 mM NaCl (pH 8.0) or 100 mM sodium phosphate, 150 mM 
        NaCl, 1 mM EDTA (pH 7.5)"

      - Oligo concentration: 1 pmol/µL (1 µM)

   - It seems like it doesn't matter too much.  I'd rather use a thermocycler, 
     because that seems much more reproducible.  The 45 minute protocol is 
     probably fine.

- "Ligation was performed by the addition of 3 U T4 Kinase and 20 U of T4 RNA 
  ligase at 25°C for 10, 20, and 40 min."

   - The Results section clarifies that they tried 10, 20, and 40 min, and 
     found that the reaction was complete after 10 min.  So I should just 
     incubate for 10 min.

   
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

Ligate linker-N --- 2019/12/9
-----------------------------
.. protocol::

   See binder: 2019/12/9 and 2019/12/13

.. figure:: 20191213_ligate_linker_n.svg

- The transcribed RNA is not very homogeneous.  I can think of a few likely 
  reasons for this:

   - Too much RNase (leftover from miniprep) in template DNA.  I could test for 
     this by doing a side-by-side comparison between the template I used this 
     time and template that I purified by phenol-chloroform extraction.

   - Expired HiScribe kit.  I remember getting smeary RNA with older IVT kits 
     in grad school, so it wouldn't surprise me if that was happening here.  If 
     none of my other ideas explain the heterogeneity, I should just try 
     ordering a new kit.

   - I contaminated something.  I could test for this by repeating the 
     transcription and just being more careful.  I thought I was pretty careful 
     this time, though...

- The ligation was 66% efficient, less than the 90--95% efficiency reported by 
  [Naimudden2016]_.  But I have a number of things I can try (discussed in the 
  Ligation_ section above) to improve this.

  Note that this efficiency is probably a slight overestimate.  I calculated 
  efficiency using the same equation as [Naimudden2016]_, but this equation 
  doesn't account for the fact that the conjugate has 28 bp of double-stranded 
  DNA/RNA hybrid.  `According to Biotium 
  <https://biotium.com/faqs/gelred-gelgreen-ssdna-rna/>`, "titration assays 
  using a fluorescence microplate reader showed that the fluorescence signal of 
  GelRed® bound to ssDNA and RNA is about half that of GelRed® bound to dsDNA."  
  Assuming that double-stranded DNA/RNA is as bright as dsDNA, this would give 
  a corrected efficiency of 64%.

  There are also reasons why this efficiency could be just plain inaccurate.  
  One is that the smeary RNA made subtracting the background rather subjective.  
  Hopefully I can improve this by getting cleaner RNA.  Another is that there 
  could be some FITC signal in the red channel.  To check for this, I need to 
  measure both the red and green channels before adding GelRed, which I didn't 
  do this time.  Note that the efficiency looks much lower in the 300 nm GelRed 
  image.  This image shouldn't have any signal from FITC (another thing I 
  should test), but it does have a smear that could be making the lower band 
  seem brighter.

- Next time I do this experiment, I should setup control reactions without 
  linker and mRNA.  This way, all three lanes would have the same amount of 
  material, which would make the gel easier to interpret.

- Linker-N runs about with the dye front.  So don't run the dye front off the 
  gel next time.  That said, I'm mostly interested in the difference between 
  the two mRNA bands, and running the gel longer might help resolve them 
  better.

- Note sure what that high-MW linker-N band is.  (It's more easily seen in the 
  "intensity level 3" image that I didn't include here.)  But it might be a 
  consequence of the lane being severely overloaded.

- I think the green scratch is caused by the EZdoc UV tray.  The laser scanner 
  images without the scratch (not shown here) were taken before I'd added 
  GelRed or imaged with the EZdoc, and the image with the scratch was taken 
  after.  I thought the scratch could also be due to something on the bottom of 
  the tip-box scratching the gel during shaking, but the scratch (vertically 
  all the way from top to bottom, rather than circular) is not really 
  consistent with that.  

Ligate Linker-N
---------------
- Image the gel using the 488 nm and 520 nm lasers (Sapphire) and the 300 nm 
  illuminator (EZdoc) before adding any GelRed.  This will allow me to be more 
  confident about overlapping signal between fluorescent channels.
