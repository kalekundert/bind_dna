*********************
Confirm PCV2 activity
*********************

Before I can use my PCV2-Zif268 fusion to test the exonuclease binding assay, I 
need to confirm that the PCV2 domain is functional.

First try
=========
.. protocol:: 20210720_confirm_pcv2.pdf 20210720_confirm_pcv2.txt

.. figure:: 20210720_confirm_pcv2.svg

Observations:

- I see no evidence of any PCV2 activity.

- Considering my results from :expt:`123`, this probably just means that I 
  didn't add enough protein.  I loaded about 15 pmol of protein in each lane of 
  this experiment, similar to the most dilute lanes in :expt:`123`.  It's 
  interesting that ≈15 pmol of f12 is easily visible with GelGreen, while ≈15 
  pmol of Cy5-labeled DNA is quite faint.

- f12 isn't perfectly clean, but it's pretty good.

Second try
==========
After remeasuring the concentration of PCV2-Zif268 via PAGE, I'm repeating this 
experiment with (hopefully) the correct amount of protein.

.. protocol:: 20210727_confirm_pcv2.txt

.. figure:: 20210727_confirm_pcv2.svg

Observations:

- Only a small amount of the DNA reacted.  

- The shifted band in this reaction similar to the shifted band that I 
  attributed to the truncated protein, e.g. :expt:`68`.  That would make sense, 
  because the PCV2 domain and the PCV2-Zif268 fusion are about the same size.

- f12 has no spacer between the ori sequence and the spacer modification.  I 
  saw in :expt:`72` that the A15 spacer (f107) works well.  The f12 spacer also 
  has a Zif268 binding site, which could possibly be interfering with the 
  reaction (although I will need it for the final binding assay).

- My protein solution froze.  According :download:`to this table 
  <prop_aq.pdf>`, a 40% (w/w) glycerol solution should freeze at -15.5°C.  I 
  don't know exactly how cold my freezer is, but I should probably go to at 
  least 50% glycerol.

  `This website`__ has more freezing point information.  It shows that I 
  probably need at least 60% (w/w) glycerol, and that 66.7% (w/w) gives the 
  lowest possible freezing point.

- Should I have used an excess of protein?  I can't remember what the original 
  reference did.

- The SYPRO Ruby staining really didn't work well.  The ladder is barely 
  visible, and I don't see any protein at all.  My fix/destain solutions aren't 
  quite right, but I don't really think that's the problem.  I've had better 
  luck with SYPRO orange; I might just stick with that going forward.

__ https://www.engineeringtoolbox.com/glycerine-boiling-freezing-points-d_1590.html


Conclusions:

- Could I find another amplicon, ideally with Cy5?
