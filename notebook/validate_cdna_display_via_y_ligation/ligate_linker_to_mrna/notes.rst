*********************
Ligate linker to mRNA
*********************

.. toctree::
   :hidden:

   linker_n_via_naimudden2016/notes
   linker_n_via_kitamura2002/notes
   linker_n_via_optimized_conditions/notes
   pseudo_linker_via_5_phosphorylation/notes
   poly_a_linker_via_optimized_conditions/notes
   poly_a_linker_via_crowding_agents/notes
   poly_a_linker_via_incubation_time/notes
   poly_a_linker_via_linker_concentration/notes


[Naimudden2016]_ reports a 90--95% efficiency ligating linker-N to mRNA.  My 
goal is to achieve similar efficiencies.

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
   which would've made the linker-N band artificially dim.  This step is not 
   mentioned in [Naimudden2016]_, but is mentioned in [Naimudden2011]_.

.. _validate_cdna_display_ligation:

Considerations
==============

Ambiguities in [Naimudden2016]_
-------------------------------
.. update:: 2020/03/31

   The ambiguities discussed below are completely resolved by [Kubo2020]_.
   
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

Other protocols
---------------
After comparing a number of published T4 RNA ligase protocols (see below), I 
think there are several things to think about trying:

- The 10 min incubation recommended by [Naimudden2016]_ is an outlier, although 
  I understand the desire to be gentle with RNA.  Still, I should try longer 
  incubation times.

  This is from the NEB FAQ for `T4 RNA ligase I 
  <https://international.neb.com/faqs/2018/01/30/what-is-the-optimal-reaction-temperature-and-time-for-t4-rna-ligase-i>`__:

    In general, increased reaction time and lowered reaction temperature yield 
    more complete ligation reactions. A typical ligation reaction should be 
    carried out at 25°C for up to 2 hours. For longer oligos, 2 hours of 
    incubation at 25°C followed by overnight incubation at 16°C may improve 
    yield.

- The protocols really differ in the amount of enzyme used:

  - KBK: 4 U/pmol
  - [Naimudden2016]_: 0.4 U/pmol
  - [Kubo2020]_: Units not specified
  - Takara: 2.5–5.0 U/pmol
  - NEB: 0.5 U/pmol
  - [Kitamura2002]_: 5 U/pmol

  That said, I'm already on the high end, so I don't think that adding more 
  ligase is likely to help much.

- All of the non-cDNA-display protocols have 25% PEG, which is known to 
  dramatically improve ligase activity [Bauer2017]_.

- NEB recommends using a 2–10x excess of linker.  That's in line with what I 
  was thinking about trying, anyways.

- Most of the protocols call for adding ATP separately from the buffer.  It's 
  possible that this is important to keep the ATP active, but I do a good job 
  taking care of my T4 DNA ligase buffer, so I don't think I need to worry 
  about this.

:expt:`17`
~~~~~~~~~~
- Setup the ligation reaction:

  - 5 pmol mRNA
  - 5 pmol linker
  - 0.1x PBS (left over from annealing reaction)
  - 1x T4 DNA ligase buffer
    - 50 mM Tris-HCl, pH 7.5
    - 10 mM MgCl₂
    - 10 mM DTT
    - 1 mM ATP
  - 0.01% BSA
  - 20 U T4 RNA ligase
  - water to 40 µL

- Incubate at 25°C for 10 min, then 65°C for 10 min

I'm not sure where exactly this protocol came from.  Some thoughts:

- 10x diluted PBS: my own optimizations of the annealing reaction.

- 20 U T4 RNA ligase: Probably [Naimudden2016]_, even though I'm using 10x less 
  mRNA/linker (I probably assumed that more enzyme wouldn't hurt, and this is 
  already about as little as I can pipet).

- BSA: Probably because Takara included BSA with the enzyme, although I'm using 
  a different concentration that their online protocol suggests.

- 10 min incubation time at 25°C: [Naimudden2016]_.  I don't know where the 
  65°C incubation/heat denaturation came from, though.

[Naimudden2016]_
~~~~~~~~~~~~~~~~
- Setup the ligation reaction:

  - 50 pmol mRNA
  - 50 pmol linker-N
  - 3 U T4 PNK (NEB)
  - 20 U T4 RNA ligase (Takara)
  - Unspecified total volume, but assuming 50 µL based on [Kubo2020]_.
  - Unspecified buffer (possibly no buffer)

- Incubate at 25°C for 10, 20, or 40 min.

[Kubo2020]_
~~~~~~~~~~~
- Setup ligation reaction:

  - 40 pmol mRNA
  - 40 pmol linker
  - 1x T4 DNA ligase buffer
  - 2 µL T4 PNK (unspecified concentration)
  - 2 µL T4 RNA ligase (unspecified concentration)
  - water to 40 µL

- Incubate at 25°C for 15 min (or at 16°C for 2h)

Takara --- T4 RNA Ligase (2050A)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`Original protocol 
<https://www.takarabio.com/assets/documents/User%20Manual/2050A_DS.v1901Da.pdf>`__.

- Setup the ligation reaction:

  - 1-2 µg ssRNA (8-16 pmol f11)
  - 1x T4 RNA Ligase buffer (same as NEB T4 DNA ligase buffer)
    - 50 mM Tris-HCl, pH 7.5
    - 10 mM MgCl₂
    - 10 mM DTT
    - 1 mM ATP
  - 0.006% BSA
  - 25% PEG 6000
  - 40-50 U T4 RNA ligase
  - water to 50 µL

- Incubate at 5-16°C for 16-18h.

- Stop the reaction by adding 2 µL 500 mM EDTA.

NEB --- T4 RNA Ligase I (M0204)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
`Original protocol 
<https://www.neb.com/protocols/2018/10/17/protocol=ligation=of=an=oligo=to=the=3=end=of=rna=using=t4=rna=ligase=1m0204>`__.  
I'm using T4 RNA ligase from Takara, but this protocol may still be relevant.

- Setup the ligation reaction:

  - 20 pmol RNA
  - 40-200 pmol DNA or RNA oligo
  - 1x T4 RNA Ligase Reaction Buffer (T4 DNA ligase buffer w/o ATP)
    - 50 mM Tris-HCl, pH 7.5
    - 10 mM MgCl₂
    - 1 mM DTT
  - 1 mM ATP (this makes the buffer equivalent to T4 DNA ligase buffer)
  - 10% DMSO (optional)
  - 15-25% PEG 8000
  - 20 U Murine RNase inhibitor (M0314, optional)
  - 10 U T4 RNA Ligase 1
  - water to 20 µL

- Incubate at 25°C for 2 hours or at 16°C for 16 hours

- Stop the reaction with a spin-column cleanup.

[Kitamura2002]_
~~~~~~~~~~~~~~~
- Setup the ligation reaction:

  - 10 pmol 5' end
  - 10 pmol 3' end
  - 50 mM Tris, pH 8.0
  - 10 mM MgCl₂
  - 0.1 mM ATP
  - 0.001% BSA (10 mg/L)
  - 1 mM hexamminecobalt (III) chloride
  - 25% PEG 6000
  - 50 U T4 RNA ligase
  - water to 10 µL

- Incubate at 25°C for 16h.
