**************************
FLAG, o237, unreacted mRNA
**************************

Previously, I've made the following observations:

- :expt:`1`: increasing the linker concentration has no effect on yield.  
- :expt:`3`: increasing reaction time has no effect on yield

Together, these observations make me wonder if there's any difference between 
the mRNAs that do and don't react.  For example, maybe the unreacted mRNAs  
have untemplated 3' additions that interfere with Y-ligation.  To test this, I 
purified the unreacted mRNA by PAGE purification (:expt:`108`) and will see 
what happens when I use this mRNA for the ligation reaction.

Results
=======

.. protocol:: 20210610_ligate_unreacted_mrna.pdf 20210610_ligate_unreacted_mrna.txt

.. figure:: 20210610_ligate_unreacted_mrna.svg

Observations:

- The unreacted mRNA from the previous ligation reaction is ligated much less 
  efficiently than the fresh mRNA.  This supports the idea that there's 
  something about the mRNA itself that is preventing the ligation.  My 
  hypothesis remains that there is untemplated addition to the 3' end.

- The bands are very faint, so don't read too much into the densiometry 
  numbers.  But the difference between the fresh and reused mRNA is 
  qualitatively clear.  I realized that I wasn't loading enough mRNA, so I 
  updated my config files to address that in the future.

- I don't know why the Cy5 bands are so dim.  I double-checked that the linker 
  should've been 1 µM in every sample, so I can't explain why the ligation 
  reactions would have so much less Cy5 than the linker-only control.
