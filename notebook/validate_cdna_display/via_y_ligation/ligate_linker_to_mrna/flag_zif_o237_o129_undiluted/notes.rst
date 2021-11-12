******************************
FLAG/Zif, o237/o129, undiluted
******************************

Previously, I found that salt is required for annealing and inhibitory for 
ligation (:expt:`12`).  However, there are reasons to revisit this conclusion:

- [Reyes2021]_ does not include salt, and I'm trying to follow this protocol.

- My methods have improved.  I now:

  - Label the linker with Cy5, to avoid cross-talk with the mRNA signal.

  - Include o194 in the electrophoresis sample buffer, to clearly distinguish 
    annealed from ligated product.

- Linker-N required phosphorylation, while my current linkers do not.  My 
  concerns was about salt interfering with phosphorylation, not ligation 
  (:expt:`59`).

If salt is not inhibitory to ligation, then I can remove the 10x dilution from 
my protocol, which would make it easier to gel purify the ligated product.

Results
=======
.. protocol:: 20210514_ligate_undiluted_salt.txt

.. figure:: 20210514_ligate_undiluted_salt.svg

Observations:

- Diluting the reaction does slightly improve ligation efficiency.

- The ligation reactions were less efficient (20-30%) than they usually are 
  (40-50%).  I don't know why; as far as I know I set up the reaction in 
  exactly the same way as I did for :expt:`111`.

- The second image I took (to reduce saturation) has a weird phenomenon in the 
  SYBR Green II channel where the centers of the bands are dimmer than the 
  edges.  This seems like some sort of artifact, so take the Zif268 densiometry 
  results with a grain of salt.

- The FLAG linker (o237) doesn't need o194.  
  
  o194 is specific to the Y-tag sequence, which f11/o237 don't have, so there 
  was really no reason to add o194 to these lanes in the first place.  Despite 
  that, there is no "ligated" band in the −ligase control.  This is in contrast 
  to the Zif268 −ligase control, which actually does have a very faint 
  "ligated" band.  My hypothesis is that the GSGS linker denatures more easily 
  because it has a lower GC content.  In any case, it's nice to know that I 
  don't need o194 for these reactions.

- The FLAG mRNA band is faint even with SYBR Green II staining.  I'll probably 
  go back to using GelGreen in the future.

- The rightmost lane is loading-dye only.  This shows that the faint band 
  around 53 nt in the SYBR Green II channel is due to the loading dye, and most 
  likely due to o194.

