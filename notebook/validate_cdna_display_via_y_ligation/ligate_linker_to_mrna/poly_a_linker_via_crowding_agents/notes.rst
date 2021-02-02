*********************************
Poly-A linker via crowding agents
*********************************

.. figure:: 20200817_optimize_peg.svg

.. datatable:: 20200817_optimize_peg.xlsx

  Densiometry of the stained gel.

.. datatable:: 20200817_optimize_peg_cy5_unsat.xlsx

  Densiometry of the unstained gel.

- Something was wrong with the GelGreen (not shown, but see other gels in data 
  directory): it wasn't very fluorescent in the expected channel, and it was 
  fluorescent in the Cy5 channel.  I made a fresh solution of 3x GelGreen for 
  this experiment, and I'm quite confident I actually used GelGreen (i.e. not 
  GelRed, which I also have).  I remember checking the label, and the solution 
  was the orange-ish color of GelGreen (not the red-ish color of GelRed).  I 
  didn't use NaCl in my 3x solution this time, but that hasn't caused problems 
  previously.

  Because the stain behaved unexpectedly, I have to interpret the bands in the 
  stained gel (left in the above figure) with caution.  Since bands for the 
  linker, the mRNA, and the fusion product are all present, I assumed that the 
  channel integrates both the mRNA and linker signal, thus I can calculate 
  yield as:
  
  .. math::
  
    \frac{{px}_{mRNA + linker}}{{px}_{mRNA} + {px}_{linker} + {px}_{mRNA + linker}}

  The caveat is that I don't really know what's being measured here, so I need 
  to take the results with a grain of salt.  That said, the densiometry results 
  are consistent with what I've seen in previous experiments.
  
  Due to the uncertainly about the GelGreen, I ran another gel using the left 
  over reaction (they were 10 µL reactions, and I only used 4.44 µL per gel).  
  Note that these reactions sat for several hours at 4°C.  This second gel I 
  did not stain, and imaged only the Cy5 channel.  This alone should be enough 
  information to determine how far the reaction progressed, as I've previously 
  seen that both the Cy5 (linker) and GelGreen (mRNA) channels give consistent 
  results.  Still, I'm more comfortable having both channels.

- The addition of PEG does not seem to have a measurable effect on the yield of 
  the reaction.  This surprises me, given that PEG is recommended in most 
  protocols and is widely known to improve the efficiency of T4 DNA ligase.  I 
  think these results are trustworthy, though.  I still might keep adding PEG 
  just because it doesn't seem to hurt.

- Both reactions appeared to have progress significantly further (75% vs 50%) 
  in the second gel, i.e. after incubating at 4°C for several hours.  This 
  suggests that a longer incubation time, even at very low temperatures, is 
  likely to be beneficial.

- 
