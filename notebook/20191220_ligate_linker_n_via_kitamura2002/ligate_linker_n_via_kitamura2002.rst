***********************************
Ligate linker-N via [Kitamura2002]_
***********************************
[Kitamura2002]_ has a better description of the Y-ligation method than any 
other paper I've found (including [Nishigaki1998]_, which introduces the idea).  
This version of the protocol is for ligating two DNA oligos (although T4 RNA 
ligase is still used), but I think some of the ideas are worth trying:

- "The 5′-half and the 3′-half DNAs in 10 μl of ligation buffer containing 50 
  mM Tris–HCl (pH 8.0), 10 mM MgCl2, 10 mg/l BSA, 1 mM hexamminecobalt (III) 
  chloride, and 25% polyethylene glycol 6000 (Tessier et al., 1986) were 
  annealed through heating at 94°C (5 min) and then at 60°C (15 min). Ligation 
  was then carried out with 50 U T4 RNA ligase (Takara) in the presence of 0.1 
  mM ATP at 25°C for 16 h."

  This ligation buffer is similar to the Takara T4 RNA ligase buffer, with the 
  following differences:
   
   - No DTT
   - 25% PEG
   - 1 mM hexaamminecobalt (III) chloride

  This buffer is fairly different from what I'm using.  According to `Wikipedia 
  <https://en.wikipedia.org/wiki/Hexamminecobalt(III)_chloride#Uses>`_, 
  "[Co(NH3)6]3+ is an unusual example of a water-soluble trivalent metal 
  complex and is of utility for charge-shielding applications such as the 
  stabilization of highly negatively charged complexes, such as interactions 
  with and between nucleic acids."  PEG is known to accelerate ligation,  
  although I wonder if it will accelerate intermolecular ligation more than 
  intramolecular ligation.

  16h at room temperature is an extremely long incubation for RNA.

- "The 5′-end of the 3′-half strand should be phosphorylated for a ligation 
  reaction."

  It seems like they're using phosphorylated primers.  I should've gotten 
  linker-N phosphorylated!  And with Cy5!  (Note, if linker-N had Cy5, I could 
  try visualizing my protein with LUCY-506, which has basically the same 
  spectrum as FITC.)

- "Equal amounts of the 5′-half and the 3′-half strands (usually 10 pmol each 
  in 10 μL) were combined and hybridized through their stem regions."

  This is pretty close to what I'm doing (5 pmol in 10 µL).


Results
=======

Cobalt and PEG --- 2020/01/23
-----------------------------
I want to start by following the [Kitamura2002]_ protocol pretty closely.  In 
addition, I setup reactions with/without cobalt and PEG, to see if either 
reagent had a discernible effect on its own.

To eliminate poor phosphorylation as a cause for poor efficiency, I used a 
phosphorylated oligo and skipped the phosphorylation step.  I don't have 
phosphorylated linker-N, so I ordered a new oligo from IDT (o93) similar to 
linker-N but lacking the puromycin arm.

.. protocol::

   See binder.

.. figure:: 20200124_ligate_o93_cobalt.svg

- I don't know what happened to the 16h gel.

- No detectable 

- It's odd that the red channel is so smeary.  That indicates that the linker 
  is getting ligated to something, but not the full length mRNA.  And there are 
  no red bands above the control f11 band, so it's not that o93 is polymerizing 
  or anything like that.

  The linker is being shifted, but the mRNA isn't really running any slower.

  .. todo::

      What if I make a smaller "mRNA" so I can see the ligation shift more 
      easily?  f11 is 388 nt.  o93 is 24 nt (noticeably smaller than linker-N).  
      Maybe something on the order of 50 bp would work well.  If I used the His 
      tag as a reverse primer to delete everything up to the T7 promoter, 
      that'd give me something about 65 nt long.  That might be worth making.

  .. todo::

      - Add RNase inhibitor.
      - Order nuclease-free reagents (e.g. NaCl, water)
      - Check for trace metal levels in cobalt salt.
      - Make sure to do 70°C step.
      - Formamide dyes generally result in better denaturing, might use that.

      - −enzyme control to make sure denaturing is working (i.e. unligated RNA 
        and DNA being separated.)
