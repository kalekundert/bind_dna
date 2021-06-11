*****************************
Visualize FLAG via FluoroTect
*****************************
I realized that the FLAG peptide is small enough that it's kinda hard to see on 
a gel.  In this experiment, I tried two things to make it easier to see:

- Include vs. exclude the tRNA digestion step

- Bis-tris vs. tricine SDS PAGE

  - Tricine gels are known to give better resolution for peptides.

Alternatively, I could try visualizing FLAG via Western blot (instead of 
FluoroTect GreenLys).  This would be very specific, but it might not be easy to 
see if the DNA and the protein are superimposed.  I'm not pursing this idea at 
the moment, but it's something to keep in mind.

Results
=======
.. protocol:: 20210415_optimize_flag_visualization.txt

.. figure:: 20210415_optimize_flag_visualization.svg

Observations:

- If looking at the raw image for the tricine gels, the notched gel is the 
  10-20% gradient.

- I can't see the 17 kDa ladder band on the tricine gels.  I don't know why.

  .. update:: 2021/05/06

    After looking at some more of these gels (e.g. :expt:`99`), I think the 14 
    and 17 kDa ladder bands are just nearly superimposed.

- The crystal violet loading dye moves into the gel, despite its positive 
  charge.  This is the case even with (i) only loading buffer and (ii) only 
  glycerol + crystal violet.  These lanes are not in the above figure, but are 
  visible in the raw data for the bis-tris gel.  I think this is probably due 
  to the SDS binding the dye.

  That said, crystal violet does seem to have less of a shadow that SERVA blue 
  G250/phenol red.  It is slightly fluorescent in the Coomassie channel, but 
  not at all fluorescent in the BODIPY/GreenLys channel.  It also doesn't seem 
  to quench green fluorescence like SERVA blue G250/phenol red does.

- FluoroTect GreenLys has a broad, low-MW smear.  I'm not sure what this 
  exactly is.
  
  - It's unaffected by RNase treatment.

  - It's somewhat fainter in the −RNase and +mRNA lanes.
  
  - It's narrower on the tricine gels than on the bis-tris gels.  The FLAG 
    peptide unfortunately falls within the smear in both kinds of gels, but 
    Zif268 would be easily outside the smear in the tricine gels.

- According to Thermo, the resolution ranges for the gels tested in this 
  experiment are roughly:

  - 16% tricine: 0-6 kDa
  - 10-20% tricine: 2-200 kDa
  - 4-12% bis-tris/MES: 3-250 kDa

- The +mRNA, −RNase condition has a faint band with slightly higher MW than the 
  tRNA.  I wonder if this is tRNA that's still attached to the peptide?

Discussion
==========
- The Bolt gels actually look like the best option in this experiment.  But I'm 
  hesitant because they didn't look nearly so good in :expt:`99` (Apr 7, 2021).  
  Maybe that could be attributed to the loading dye, though.

- The RNase treatment doesn't help with visualizing FLAG.  It doesn't eliminate 
  the smear, and in fact makes it a little brighter, which makes the FLAG band 
  harder to see.  RNase treatment would be helpful if my product were just 
  slightly bigger, though.

- I don't think there's a clear reason to pick any of these gels over the 
  others.  I might revisit this once I start trying to attach mRNA.

