***********************************
Ligate linker-N via [Kitamura2002]_
***********************************
[Kitamura2002]_ has a better description of the Y-ligation method than any 
other paper I've found (including [Nishigaki1998]_, which introduces the idea).  
There are a number of differences between this protocol and what I've tried in 
:expt:`20191219_anneal_linker_n_with_salt`:

- "The 5′-half and the 3′-half DNAs in 10 μl of ligation buffer containing 50 
  mM Tris–HCl (pH 8.0), 10 mM MgCl2, 10 mg/l BSA, 1 mM hexamminecobalt (III) 
  chloride, and 25% polyethylene glycol 6000 (Tessier et al., 1986) were 
  annealed through heating at 94°C (5 min) and then at 60°C (15 min). Ligation 
  was then carried out with 50 U T4 RNA ligase (Takara) in the presence of 0.1 
  mM ATP at 25°C for 16 h."

  This buffer is fairly different from what I'm using.  According to `Wikipedia 
  <https://en.wikipedia.org/wiki/Hexamminecobalt(III)_chloride#Uses>`_, 
  "[Co(NH3)6]3+ is an unusual example of a water-soluble trivalent metal 
  complex and is of utility for charge-shielding applications such as the 
  stabilization of highly negatively charged complexes, such as interactions 
  with and between nucleic acids."  PEG is known to accelerate ligation,  
  although I wonder if it will accelerate intermolecular ligation more than 
  intramolecular ligation.

- "The 5′-end of the 3′-half strand should be phosphorylated for a ligation 
  reaction."

  It seems like they're using phosphorylated primers.  I should've gotten 
  linker-N phosphorylated!  And with Cy5!  (Note, if linker-N had Cy5, I could 
  try visualizing my protein with LUCY-506, which has basically the same 
  spectrum as FITC.)

- "Equal amounts of the 5′-half and the 3′-half strands (usually 10 pmol each 
  in 10 μL) were combined and hybridized through their stem regions."

  This is pretty close to what I'm doing (5 pmol in 10 µL).


