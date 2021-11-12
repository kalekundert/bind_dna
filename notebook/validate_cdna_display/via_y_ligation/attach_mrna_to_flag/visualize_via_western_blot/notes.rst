**************************
Visualize via Western blot
**************************

[Reyes2021]_ used western blotting to visualize the FLAG peptide.  This 
approach has some advantages and disadvantages relative to FluoroTect:

Advantages:

- Only visualize the protein, no tRNA band to ignore.
- Compatible with green DNA labels/stains.

Disadvantages:

- Slower.
- Need to composite manually.

Considerations
==============

Chemiluminescence vs. fluorescence
----------------------------------
Advantages of chemiluminescence:

- Signal amplification, may be more sensitive

  - According this this article, fluorescence (NIR, not visible) actually has 
    greater sensitivity:

    https://www.cytivalifesciences.com/en/us/Solutions/Protein-Research/Knowledge-center/western-blotting/Fluorescence-vs-Chemiluminescence

- Fluorescence

  - Quantitative, can compute fraction reacted FLAG

I think NIR fluorescnce is the way to go (although [Reyes2021]_ did 
chemiluminescence).

Primary antibody
----------------
I'm just going to use the same anitbody as [Reyes2021]_; presumably that will 
at least be a reasonable choice:

- Anti DYKDDDDK tag, Monoclonal Antibody (FUFIFILM Wako Pure Chemical, 
  018-22381)

[Reyes2021]_ specifically used a HRP-fused antibody (i.e. no secondary), but 
fortunately FUJIFILM seems to sell an unlabeled version of this antibody.

Secondary antibody
------------------
The primary described above is mouse, so I need anti-mouse.

Use IRDye 800CW for detecting the protein of interest, and IRDye680RD for 
detecting housekeeping genes (which I don't have).

- https://www.abcam.com/secondary-antibodies/irdye-conjugated-secondary-antibodies
- The specific secondary from this example is goat anti-mouse (Abcam ab216772), 
  so I think it would work for me:

  https://www.abcam.com/goat-mouse-igg-hl-irdyereg-800cw-preadsorbed-ab216772.html

Membrane
--------
Pore size:

- 0.2 µm and 0.45 µm are the typical options.

- 0.45 µm is the "generally suitable" option.

- 0.2 µm membrane are recommended for proteins with MW < 20 kDa, which 
  definitely includes FLAG.

  https://www.licor.com/bio/guide/westerns/which_membrane

Polymer:

- [Reyes2021]_ used PVDF.

- PVDF and nitrocellulose both bind proteins via hydrophobic interactions.

  https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/western-blot-transfer-methods.html

- Nylon membranes bind by electrostatic/ionic interactions.

- Since FLAG is very polar, nylon might actually work the best.

- I'll start with PVDF in the interest of following [Reyes2021]_.  But I also 
  have nylon membranes on hand, so I might be able to try them as well.

Blocking agent
--------------

- [Reyes2021]_ used "PVDF blocking reagent" for "Can Get Signal" (which is a 
  chemiluminescence visualization kit).

- Seems like I can use either non-fat milk or BSA, both of which I have on 
  hand.

  https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/blocking-buffers-western-blot-elisa.html

- Can use a blocking buffer specifically designed to limit background 
  fluoresence, although milk seems to have low background on its own.  (BSA 
  presumably also has low background.):

  `Blocker™ FL Fluorescent Blocking Buffer, 10X (Thermo 
  37565) <https://www.thermofisher.com/order/catalog/product/37565#/37565>`_

- I think I'll start with milk.  It's cheap, and should be fine.

iBlot vs. semidry transfer
--------------------------
The Wyss has an iBlot1 machine (really it's just called iBlot, but I'll refer 
to it as iBlot1 to distinguish it from iBlot2).  Dima has two iBlot2 machines, 
but they're specifically for his project (I asked him).  So if I decide to do a 
dry transfer, I'd want to get iBlot1 stacks.

iBlot1 PVDF 0.2 µM mini-gel stack:

- Catalog #: IB401002
- Price: $225

Advantages of iBlot:

- Fast and easy
- Comparable in price to buying all reagents for semi-dry

Advantages of semi-dry:

- Would be cheaper in bulk
- What [Reyes2021]_ did.
- I think this is considered higher quality, but not sure about that.
- Can optimize buffers.

I'm going to go with iBlot:

- It's the easiest to get started with.
- I don't plan to be doing lots of gels.
- I think the tranfer will work well; I won't need to be optimizing buffers 
  etc.

Transfer buffer
---------------
https://www.bio-rad.com/en-us/applications-technologies/types-western-blot-transfer-buffers?ID=LUSQA88UU

Buffers specifically designed for semidry transfers have higher ionic strength, 
the higher voltage/current during the transfer process (relative to wet 
transfer) will more quickly exhaust the buffer.  That said, it seems like 
people can still get good results using "regular" buffers for semidry 
transfers.

Buffer components:

- Buffer, pH usually seems to be quite basic, e.g. 9-11.

- SDS:

  "SDS and alcohol play opposing roles in a transfer. SDS in the gel and in the 
  SDS-protein complexes promotes elution of the protein from the gel but 
  inhibits binding of the protein to membranes. In cases where certain proteins 
  are difficult to elute from the gel, SDS may be added to the transfer buffer 
  to improve transfer. SDS in the transfer buffer decreases the binding 
  efficiency of protein to nitrocellulose membrane; PVDF membrane can be 
  substituted for nitrocellulose when SDS is used in the transfer buffer. 
  Addition of SDS increases the relative current, power, and heating during 
  transfer and may affect the antigenicity of some proteins."

- Alcohol: 

  "Alcohol (methanol or ethanol), on the other hand, removes the SDS from 
  SDS-protein complexes and improves the binding of protein to nitrocellulose 
  membrane but has some negative effects on the gel itself. Alcohol may cause a 
  reduction in pore size, precipitation of some protein, and some basic 
  proteins to become positively charged or neutral. All of these factors will 
  affect blotting efficiency."

  "Only high-quality, analytical grade methanol should be used in transfer 
  buffer; impure methanol can increase transfer buffer conductivity and result 
  in poor transfer."

