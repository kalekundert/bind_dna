*********
First try
*********

Considerations
==============

Controls
--------
- An ideal negative control for testing the puromycin coupling reaction would 
  be o129 without puromycin.  But that would cost $400, and would really only 
  be useful for this experiment.  I might not even need it (i.e. if the 
  reaction clearly works), so I'll proceed without for now, and think more 
  about ordering this oligo (o197) if I need it.

- In :expt:`19` (2020/02/18), my "controls" were to leave out reactions.  This 
  is nice because it also allows me to monitor the reaction and see where 
  problems arise (if there are problems).

  I do want to include a "−ligation, +expression" sample (which is not 
  represented in the aforementioned experiment, because I think it might be 
  informative to see mWasabi covalently coupled to just the puromycin linker, 
  but not the mRNA.

Volumes
-------
I'm trying to work out how much mRNA to load on the gel.  This is both to 
decide how much o194 to add the loading buffer, and the work out how big to 
make each reaction.  I'm using :expt:`19` as a reference:

- Assuming that I concentrated the mRNA to the maximum ability of the Amicon 
  spin filter (15 µL), I added 333 nM mRNA to the expression reaction.  Note 
  the ideal mRNA concentration as indicated by :expt:`18` is 1 µM, so I'm in 
  the right ballpark.

  .. math::

    \pu{333 nM} = \frac{\pu{40 µL} \times \pu{125 nM}}{\pu{15 µL}}

- The concentration of the mRNA in the:

  - expression reaction: 26.6 nM
  - loading buffer: 12.4 nM

- The gel was overloaded, though.  In my notes, I thought it'd be better to add 
  only about half the sample, which would be 6.8 nM mRNA.

- I can do this reaction in just the same way, I just need to keep the volume 
  of the ligation reaction the same.

.. note::

   At this concentration, I don't know if I'd be able to see the mRNA even if I 
   stained for it.  But I'll keep the mRNA-only lane, because it could be 
   useful.

Ladders
-------
Since I'm not planning to stain my gel, it doesn't really make sense to run 
ladders (although I might anyways, in case I decide to post stain).  But for 
the future, I might also think about ordering a fluorescent protein standard, 
e.g. Invitrogen LC5928.  This is basically a normal ladder with AlexaFluor 
coupled to all the proteins.  It's kinda expensive ($400/125 µL), so I'll have 
to think if I'll get enough use out of it to be worth it.

Results
=======

2020/08/25:

.. protocol:: 20200825_display_mrna.txt
.. figure:: 20200825_attach_mrna_to_mwasabi.svg

- The protein expression reaction appears to be degrading the mRNA.  The 
  "annealed but unligated" and "ligated and filtered" lanes both appear to 
  contain full-length mRNA, while the corresponding lanes after protein 
  expression both do not.

  This may be a sign that I incubated the expression reaction for too long.  
  [Barendt2013]_ incubates the PURExpress reaction for 30 min at 37°C, then 10m 
  at room temperature.  I'd be a little surprised if this made a big 
  difference, but it's probably worth doing a time course.

- Note that the low-MW Cy5 band is shifted in the +everything lane as compared 
  to all of the other o129 bands.  I think this is evidence in favor of 
  nuclease contamination.  As illustrated below, the idea is that RNase H would 
  only cleave the RNA that was duplexed with DNA.  This would leave at least 4 
  un-paired RNA nucleotides (i.e. those making up one "arm" of the Y-ligation, 
  and maybe a few more 5' of that) attached to the linker, which might be 
  enough to explain the difference in apparent band sizes.

  .. figure:: band_size_hypothesis.svg

- There is no indication that any of the puromycin reacts with any of the 
  protein.  Even if the mRNA ends up being degraded (as discussed above), the 
  Cy5 would stay coupled to the peptide and be visibly retarded if the 
  puromycin reaction occurred.

  I'm pretty sure that IDT shipped me the right oligo.  I checked the IDT mass 
  spectrometry quality control data, and it seems correct:

  - Expected MW: 11236.1 Da
  - ESI peak:    11235.7 Da

  Some things I could try:

  - Incubate with high concentrations of Mg²⁺ and K⁺.  See :expt:`62`.

  - Use the Spacer-18 linker.  See :expt:`61`.

  - Run a control with the mRNA (f85) and only the puromycin arm of the linker 
    (o125, o126).  The puromycin arm---basically free puromycin---really should 
    react with the protein.  I'd expect to see something of a smear in the Cy5 
    channel, as the puromycin truncates the mWasabi gene at different points.  
    The higher-MW parts of the smear might also be fluorescent in the GFP 
    channel, if they contain the matured fluorophore.

  - Use a gene from [Barendt2013]_.  This is assuming that there's something 
    problematic about GFP, which just seems unlikely.  Also [Barendt2013]_ used 
    ankyrin repeats, which I understand can be hard to work with.

    .. note::

       Most (all?) mRNA display protocols recommend using radioactive 
       methionine for protein expression.  This would certainly make the 
       protein easy to visualize, but it would also make everything harder to 
       work with.  I think my approach of using fluorescent tags is better.  
       For now I'm using mWasabi so I can see the protein, but once I have a 
       protocol worked out, I can use non-fluorescent proteins and just monitor 
       the display reactions via the Cy5 in the linker.

- I don't know why the "annealed but unligated" lane appears to have a 
  significant amount ligated product.  I ran exactly this lane in :expt:`50` 
  and saw no mRNA/linker band at all.  Maybe it's possible that the difference 
  has to do with this being an SDS gel, rather than a urea gel, but that 
  doesn't seem particularly likely. 

- The filtration step doesn't seem particularly effective at removing unligated 
  linker, and also seems to lose a significant amount of material.  Adding 
  competing oligo (o194) might help with the first problem.  The second problem 
  may be a fluke, because I got good yield in :expt:`19`.

