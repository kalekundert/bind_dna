***********************
Detect RNase H activity
***********************

2020/11/06

George hypothesized that the degradation I'm seeing in my in vitro expression 
reactions (e.g. :expt:`66`) maybe due to RNase H contamination.  To investigate 
that possibility, I want to find an assay that I can use to test specifically 
for this activity.

Note that if I do find RNase H activity, I have a few options:

- Try to inhibit it, e.g. by adding an excess of decoy substrate or gentamycin 
  [Zhao2017]_.
- Find a kit that doesn't have any RNase H activity, e.g. maybe PUREfrex
- Design my mRNA using splint ligation, to avoid creating a RNA/DNA duplex.
- Make the dsDNA before expressing the protein.

Considerations
==============

Spectrophotometers
------------------
I want to be sure that I can actually measure the output from any of these 
assays, so I checked the wavelengths supported by the various 
spectrophotometers in the lab.

- Nanodrop 2000c: 190-840 nm

- Synergy HT: 200-999 nm

- Synergy H1: unsure, may depend on installed lasers/filters

- SpectraMax M5: 250-850 nm

Expected activity
-----------------
According to NEB, 1 unit of RNase H is the amount of enzyme necessary to 
completely digest 20 pmol of 50 bp DNA/RNA duplex in 20 min at 37°C.  My 
reactions have 1 µL of 333 nM mRNA, which is 0.3 pmol.  Since the whole 
reaction has a volume of 3.6 µL, a RNase H concentration of 4.2 U/mL would be 
sufficient to digest all of my mRNA.  I'd probably like to be able to detect 
10x less than this, just to be confident I'm mot missing anything.

Assays
------
Below is a brief review of some RNase H assays from the literature.  Having 
done this review, I'm inclined to use a fluorophore/quencher assay with some 
sort of amplification.  Fluorophore/quencher assays are more expensive 
(compared to G-quadruplex assays), but more likely to work with complex 
mixtures.  Amplification adds complexity (which I'd rather avoid), but seems 
necessary to comfortably detect the level of activity that would be sufficient 
to digest all of my mRNA.  Given this, I'm planning the try [Wang2017]_ first.

- [Rizzo2002]_:
  
  - Detection: FAM/Dabcyl
  - Amplification: None
  
  How it works:
    
  - Start with DNA/RNA hybrid molecular beacon (i.e. hairpin that has a 
    fluorophore/quencher pair on the ends).
  - RNase H can digest the hairpin, separating the fluorophore from the 
    quencher and allowing it to fluoresce.

  Sensitivity:

  - 20 U/mL (Figure 2)

  Advantages:

  - Simple, because there's only one reagent.

- [Hu2010]_:

  - Detection: G-quadruplex
  - Amplification: None
  
  How it works:
  
  - Start with a DNA/RNA duplex and N-methyl mesoporphyrin IX (NMM).
  - RNase H can digest the RNA, revealing a single-stranded DNA oligo that 
    can fold into a G-quadruplex.
  - NMM binds the G-quadruplex and becomes fluorescent.

  Sensitivity:

  - ≈1 U/mL.  Figure 3 shows traces with RNase H concentrations between 0.2–4.0 
    U/mL, but the individual traces aren't labeled.

  Advantages:

  - Cheap, because no modified nucleotides are required.

- [Wang2017]_

  - Detection: FAM/BHQ1
  - Amplification: DNAzyme

  How it works:

  - Start with DNA/RNA hybrid (encoding a sequestered DNAzyme) and DNA 
    molecular beacon.
  - RNase H can digest the RNA, releasing DNAzyme.
  - The DNAzyme can cleave many molecular beacons, leading to an increase in 
    fluorescence.

  Sensitivity:

  - 0.2 U/mL (Figure 5).  The data doesn't include the −RNase H negative 
    control, but based on Figure 1, I think it'd be hard to distinguish the 
    0.06 and 0.02 U/mL traces from that control.

  Advantages:

  - Sensitive, because the signal is amplified.
  - Tested with cell-free extracts, which is pretty much exactly what I want to 
    do.

- [Wu2018]_:
  
  - Detection: G-quadruplex/ThT
  - Amplification: Nicking restriction enzyme

  How it works:
    
  - Start with a DNA/RNA duplex, a DNA hairpin, and thioflavin-T (ThT).
  - RNase H can digest the RNA, leaving a single stranded DNA oligo.
  - This oligo can bind to the hairpin, revealing a restriction site for a 
    nicking enzyme.
  - When this site is nicked, an oligo that can fold into a G-quaruplex is 
    freed from the hairpin, and the original DNA oligo is released to bind 
    another hairpin.
  - Thioflavin-T binds the G-quadruplex and becomes fluorescent.
  - Note that the signal is amplified because each digested RNA/DNA duplex 
    can give rise to many G-quadruplexes.

  Sensitivity:

  - Can't access that part of the manuscript, but the abstract claims 0.03 
    U/mL.  The other claims I've seen have been a bit exaggerated, so I'd guess 
    this is more like 0.1 U/mL.

  Advantages:
  
  - High sensitivity, because signal is amplified.
  - Cheap, because no modified nucleotides are required.

  Concerns:

  - Can't access the whole article...
  - Thioflavin-T is also used to stain amyloid fibrils, and is not considered 
    to be a very specific binder [Wikipedia].  So this assay may not work as 
    well in a complex mixture.

- [Jung2019]_

  - Detection: TaqMan probes
  - Amplification: Taq polymerase

  How it works:

  - Start with aptamer that has a RNA/DNA stem and binds Taq polymerase.
  - RNase H can digest the RNA, freeing Taq polymerase.
  - Taq can then extend primers bound to a ssDNA template, degrading TaqMan 
    probes in the process.

  Sensitivity:

  - 1 U/mL

- [Zhang2019]_

  - Detection: Spinach/DFHBI
  - Amplification: ?

  How it works:

  - Start with RNA/DNA hybrid hairpin.

    - The "DNA" in this case is not actually DNA, but some modified backbone 
      with increased cytoplasmic stability for in vivo assays.  But I'm just 
      going to refer to it as DNA, since that's how it behaves.

  - RNase H can digest the RNA, leaving a DNA oligo
  - The DNA oligo can hybridize with a sequestered spinach aptamer via a 
    toehold.
  - RNase H can then digest the sequestered aptamer to free Spinach.
  - Spinach binding to DFHBI can be detected.

  Sensitivity:

  - 0.005 U/mL
  - Linear from 0.005-100 U/mL, which far exceeds every other method.
  - I'm actually kinda skeptical of these values, because they're just too 
    good.

  Concerns:

  - I don't know if this really counts as amplifying, because RNase H is 
    involved in the "amplification" such that there's still only 1 signal 
    molecule per RNase event.  It seems that the initial hybrid hairpin could 
    just be left out, and the whole assay could be done with the sequestered 
    spinach aptamer.

- [Zhao2017]_

  - Detection: FAM/GO
  - Amplification: None

  How it works:

  - Start with DNA-RNA-FAM oligo that forms a hairpin.
  - The hairpin binds graphene oxide (GO) via π-stacking.
  - Graphene oxide quenches FAM fluorescence.
  - RNase H can digest the RNA, freeing the RNA from the surface, and 
    allowing FAM to fluoresce.

