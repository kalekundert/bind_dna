********************
Optimize repA linker
********************

If I am in fact having problems solubly expressing repA-fusions, the reason may 
have to do with the linker connecting the two proteins.  To investigate this 
possibility, I plan to clone and test a handful of linkers with distinct 
properties (e.g. length, rigidity).

Considerations
==============

Linkers
-------
After reading some of the literature on fusion-linkers [Chen2013b]_, I came to 
feel that the most defining feature of a fusion linker is whether it's flexible 
or rigid.  With this in mind, I chose to test a set of 6 linkers:

- [Odegrip2004]_ (Pro): GGGGSAAP

  This is the linker I've been using in all of my constructs so far.  
  
- [Odegrip2004]_ (Ala): GGGGSAAA
  
  While thinking about this experiment, I realized that it has an unintended 
  A8P mutation.  The cause is a G→C mutation in the supplement of 
  [Odegrip2004]_.  This supplement also had 4 other G→C mutations in the repA 
  coding sequence.  I found (and corrected) the repA mutations while I was 
  cloning my first CIS-display plasmids, but I didn't notice the linker 
  mutation.
  
  Ultimately, I don't think an A→P mutation is a linker sequence will be much 
  of a problem (both A and P are abundant in natural linkers), but it might 
  still be worth correcting.

- [Guilinger2014]_ ("XTEN"): SGSETPGTSESATPES

  XTEN is a repetitive ≈300 aa polypeptide that is believed to have a similar 
  effect to PEGylation.  [Guilinger2014]_ found a 16-residue consensus repeat 
  within XTEN and employed it as a flexible fusion linker.  For their 
  application (fusing Cas9 to FokI), it performed well.

- [Waldo1999]_: GSAGSAAGSGEF

  I believe this sequence was designed by hand.  It performs equivalently to 
  GGGS×4, but is shorter and less repetitious.  Both features make it easier to 
  work with.

- [Arai2001]_ ("HL4"): CTGGCAGAAGCAGCGGCGAAAGAGGCGGCGGCTAAGGAAGCCGCGGCAAAGGAAGCAGCCGCTAAAGCCGCTGCA

  The HL4 linker was designed to be helical and rigid.  From [Arai2001]_:

    The linkers, which were expected to form a monomeric hydrophilic α-helix, 
    were designed according to the previous study on a short peptide forming a 
    monomeric α-helix (Marqusee and Baldwin, 1987). In that study, the best 
    helix-forming peptide, AEAAAKEAAAKEAAAKA [(i+4)E,K], which was stabilized 
    by Glu––Lys+ salt bridges, showed ~80% helicity.

  [Arai2001]_ also tested shorter and longer linkers.  I decided to use HL4 for 
  a few reasons:

  - It was one of the authors' recommendations.

  - The inter-domain distance estimated via FRET for this linker was 57.1Å, 
    which was significantly different than the baseline distance of 50.0Å for 
    an ``LAAA`` linker.  The shorter linkers were less significantly different, 
    and the longer linkers would've been more difficult to clone.
