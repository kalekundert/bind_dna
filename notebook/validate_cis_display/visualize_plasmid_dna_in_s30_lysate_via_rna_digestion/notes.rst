*****************************************************
Visualize plasmid DNA in S30 lysate via RNA digestion
*****************************************************
In my initial experiments (:expt:`24`, :expt:`25`), I've seen that S-30 
expression lysates seem to degrade linear DNA templates.  One way to address 
this problem is to express protein from plasmid templates instead, as plasmids 
are not vulnerable to exonucleases.

The primary issue with this approach is that I do not have a good way to 
visualize plasmid DNA.  Here I will experiment with different ways to do this:

- Direct staining with GelRed or GelGreen.  This will work if the plasmid is 
  sufficiently bright against the background of any other nucleic acids (e.g. 
  ribosomes, tRNAs, etc.) in the lysate.

- Digest RNA, either enzymatically or chemically, to reduce background signal.  
  This is based on the fact that most (if not all) of the background signal 
  should come from RNAs, e.g. ribosomes, tRNAs, etc.  I don't think the lysate 
  should have any DNA besides what I add.  Chemical digestion will likely 
  denature any protein or DNA in the reaction, and so would only be useful for 
  seeing whether or not the DNA survived, not whether it was bound by repA.

- Southern blotting to specifically visualize my plasmid.

This experiment will focus on the first two approaches.  See :expt:`28` for the 
third.

Considerations
==============

How much RNase?
---------------
I can't find any estimates of how much RNA is in NEBExpress or Promega S-30 
lysate.  But the PURE system (the original, not NEB's PURExpress system) has 
32.4 µg of ribosomes and 5.6 µg tRNAs, for 38.0 µg total RNA, in a 50 µL 
reaction.  I'll assume that the concentration of RNA in the S-30 lysates is 
comparable to this, e.g. 1 µg/µL.

`Millipore Sigma 
<https://www.sigmaaldrich.com/technical-documents/protocols/biology/roche/rnase-dnase-free.html>`__ 
explains how much RNase is needed to digest a given quantity of RNA.  I can 
probably make these numbers apply to other vendors as well, because everyone 
seems to use Kunitz units for RNase A.

   0.1 mU RNase, DNase-free degrades 1 μg RNA in 30 min at 37°C in a reaction 
   volume of 50 μl PCR grade water. The protein concentration of RNase, 
   DNase-free is 0.5 μg/μl. The specific activity of the enzyme is 30 U/mg, 
   corresponding to 1.5 mU/μl; one microliter of the RNase preparation is 
   sufficient to completely degrade 15 μg RNA in 30 min. Because the exact 
   concentration of RNA in biological samples is not known exactly, and one 
   typically wants to degrade the RNA quantitatively, a excess of RNase is 
   recommended.

Based on this recommendation, I'd want to use 1.0 µL of RNase A in a 10 µL 
reaction.  That's a 15x excess, which gives some leeway for my assumptions to 
be wrong.  It also keeps the volume of enzyme below 10% of the whole reaction, 
to avoid having too much glycerol.

Salt concentration
------------------
Interesting quote from `Invitrogen 
<https://www.thermofisher.com/order/catalog/product/EN0531#/EN0531>`__:

   At low salt concentrations (0 to 100 mM NaCl), RNase A cleaves 
   single-stranded and double-stranded RNA as well the RNA strand in RNA-DNA 
   hybrids. However, at NaCl concentrations of 0.3 M or higher, RNase A 
   specifically cleaves single-stranded RNA.

I don't know how much salt is in the lysates, but I should avoid adding any 
more.

EDTA concentration
------------------
[Miall1969]_ suggests that a ≈4x excess of EDTA will most effectively denature 
the ribosome.  The PURE buffer, again assuming that it is representative, has 9 
mM MgOAc and 0.5 mM CaCl₂.  That corresponds to an ideal EDTA concentration of 
38 mM.  I might round that up to 40 or 50 mM.

Temperature
-----------
[Miall1969]_ suggests that the ribosome is 50% unfolded at:

- 65°C in 0 mM EDTA
- 50°C in 10 mM EDTA.  

`Sigma 
<https://www.sigmaaldrich.com/life-science/metabolomics/enzyme-explorer/learning-center/nucleases.html>`__ 
claims that:

   The optimal temperature for [RNase A] activity is 60 °C, although the enzyme 
   does exhibit activity from 15-70 °C.

I think it'd be best to do the reaction at 60°C, with EDTA to help loosen up 
the ribosome.

GelRed vs. GelGreen
-------------------
Normally I use GelGreen for staining nucleic acids, but in this case GelRed 
might be better.  I prefer GelGreen because it is more excited by the blue 
light used by the laser scanner.  In this experiment, though, I will also be 
expressing mWasabi, which is indistinguishable from GelGreen in terms of 
fluorescence.  It may still be hard to distinguish mWasabi and GelRed, but 
GelRed should be much better excited by the UV lamp on the GelDoc (300 nm).  
GelRed also emits at higher wavelengths than mWasabi, although I don't know 
what filter (if any) the Gel Doc uses.  Likely there will still be some 
cross-talk, but with controls I might be able to see what's going on.
