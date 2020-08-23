*****************************
Optimize linker concentration
*****************************
Background: :expt:`-1`

Although the annealing and ligation reactions seem to work with the poly-A 
linker (o129), only about half of the mRNA is modified.  An obvious way to try 
modifying more of the mRNA is to use more linker.

Considerations
==============
[Naimudden2016]_ reports that other groups have used between 2–200x excess of 
linker for Y-ligation (refs 13, 20, 26).  In my own survey of :expt:`T4 RNA 
ligase protocols <-1>`, I found that most protocols recommend no excess of 
linker (i.e. 1:1 RNA:linker), but one (from NEB) recommends a 2-10x excess.

Disadvantages of using excess linker:

- Needs to be removed, because any extra puromycin will interfere with the 
  protein expression reaction.  However, the 100 kDa MWCO spin column should do 
  this pretty effectively (o129 is about 20 kDa).

- The linker is one of the most expensive reagents.

Considering the  above, I want to limit myself to a fairly narrow range of 
linker concentrations.  Specifically, I'll start by testing:

- 1x
- 2x
- 4x
- 8x
