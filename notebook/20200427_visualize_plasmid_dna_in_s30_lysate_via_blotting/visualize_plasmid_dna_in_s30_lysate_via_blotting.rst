************************************************
Visualize plasmid DNA in S30 lysate via blotting
************************************************

See :expt:`20200421_visualize_plasmid_dna_in_s30_lysate_via_rna_digestion` for 
background information.  This experiment will focus on using Southern blotting 
to visualize plasmid DNA added to S-30 lysate.

.. todo::
   - Southern blot controls
      - lysate only
      - pure plasmid
      - 0h reaction
      - 2-4h reaction

Considerations
==============

Probe dye
---------
For Southern blotting, it's important to consider the autofluorescence of the 
membrane when choosing a probe dye.  Most membranes have high background 
fluorescence for visible wavelengths, which makes it hard to use dyes such as 
FITC and Cy5.  Instead, near-IR dyes such as AlexaFluor 680, IRDye 700, and 
IRDye 800 are preferred.  The IRDye family is specifically designed for use in 
Li-Cor imagers, but should also be compatible with most laser scanners.

I'm going to use IRDye800.  It does seem as if people generally prefer to use 
AlexaFluor 680/IRDye 700, but I don't know if the lab has a Li-Cor, and I'm 
worried that the 658 nm laser that the Sapphire would use to excite these dyes 
may trigger autofluorescence (since the laser is also used to excite Cy5, and 
Cy5 is considered to have problems with autofluorescence).

If I'm able to get Southern blotting to work, I may eventually want to use the 
technique to visualize DNA and protein at the same time.  Using IRDye800 for 
the DNA would allow me to use miRFP or one of its variants (with the biliverdin 
cofactor) for the protein.
