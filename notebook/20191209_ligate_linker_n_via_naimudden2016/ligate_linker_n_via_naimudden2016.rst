************************************
Ligate linker-N via [Naimudden2016]_
************************************
[Naimudden2016]_ reports a 90--95% efficiency ligating linker-N to mRNA.  I 
want to achieve similar efficiencies.

.. note::

   I believe that the efficiency calculated by [Naimudden2016]_ is incorrect 
   and meaningless.  The issue is that they stain with SYBR green, which has 
   almost the same spectrum as the FITC in linker-N.  This makes the product 
   band disproportionately brighter in a way that cannot be accounted for.

   [Naimudden2016]_ doesn't specify what laser and filter is used to image the 
   gel (they just say "Ligation efficiency was checked by denaturing 
   polyacrylamide gel electrophoresis (6% gel and 8 M urea) using FITC and Sybr 
   Gold (Molecular Probes, USA) staining on a fluoroimager (Bio-Rad, Hercules, 
   CA, USA). Ligation efficiency was calculated by the equation: % ligation = 
   A/(A+B) × 100%".  But I don't think there's any way they could've used SYBR 
   green and not had this problem.

   It is worth noting that the left gel in Figure 1d, which shows only FITC, 
   still makes it seem like a majority of the FITC was ligated to the mRNA.  
   But they may have done an rNA spin-column clean up before running that gel, 
   which would've made the linker-N band artificially dim.  (This step is not 
   mentioned in [Naimudden2016]_, but is mentioned in [Naimudden2011]_).

.. _validate_cdna_display_ligation:

Considerations
==============
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

      - Drop dialysis after annealing to remove salt.  The half-life of NaCl 
        during drop dialysis is 15 min [Gorisch1988]_.

      - Do the annealing reaction in a small volume (e.g. 5 µL), then dilute 
        for the ligation reaction (e.g. to 50 µL) to reduce the salt to a 
        non-inhibitory level.

      - Just use T4 PNK buffer and don't worry about salt.

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
     <https://www.sigmaaldrich.com/technical-documents/protocols/biology/annealing-oligos.html>`__ 
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
     <https://www.idtdna.com/pages/education/decoded/article/annealing-oligonucleotides>`__ 
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
     <https://tools.thermofisher.com/content/sfs/brochures/TR0045-Anneal-oligos.pdf>`__ 
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
.. protocol:: 20191209_transcribe_rna.txt

   See binder: 2019/12/9 and 2019/12/13

.. figure:: 20191213_ligate_linker_n.svg

- The transcribed RNA is not very homogeneous.  See 
  :expt:`20191216_optimize_mrna_transcription` for more discussion.

- The ligation was 66% efficient, less than the 90--95% efficiency reported by 
  [Naimudden2016]_.  But I have a number of things I can try (discussed in the 
  :ref:`validate_cdna_display_ligation` section) to improve this.

  Note that this efficiency is probably a slight overestimate.  I calculated 
  efficiency using the same equation as [Naimudden2016]_, but this equation 
  doesn't account for the fact that the conjugate has 28 bp of double-stranded 
  DNA/RNA hybrid.  `According to Biotium 
  <https://biotium.com/faqs/gelred-gelgreen-ssdna-rna/>`_, "titration assays 
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
  "intensity level 3" image that I didn't include here.)  But it's also visible 
  in [Naimudden2016]_.

- I think the green scratch is caused by the EZdoc UV tray.  The laser scanner 
  images without the scratch (not shown here) were taken before I'd added 
  GelRed or imaged with the EZdoc, and the image with the scratch was taken 
  after.  I thought the scratch could also be due to something on the bottom of 
  the tip-box scratching the gel during shaking, but the scratch (vertically 
  all the way from top to bottom, rather than circular) is not really 
  consistent with that.  

