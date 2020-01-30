**********************
Validate pseudo-linker
**********************

I designed a "pseudo-linker" (o93) to use while optimizing the ligation 
reaction.  The pseudo-linker is basically Linker-N with the following 
differences:

- No puromycin arm, so I could order it from IDT for a reasonable price.
- 5' phosphorylated to (possibly) improve ligation efficiency.
- Cy5-labeled instead of FITC-labeled to facilitate 2-channel gel imaging.

Here I want to compare the pseudo-linker with linker-N, to verify that the 
former is an appropriate proxy for the latter.

Results
=======
.. protocol:

   See binder, 2020/01/28

.. figure:: 20200128_o93_vs_linker_n_9_4.svg

.. datatable:: o93_vs_linker.xlsx

Discussion
==========
- I experimented with the intensities of the lasers on the Sapphire to find 
  intensities that were as bright as possible without saturating any pixels.  
  These intensities are shown in the table above.  I should be able to continue 
  using these intensities as long as I continue using 0.5 µM linker/mRNA in 
  these reactions.

- linker-N has 17.3x less non-specific annealing than o93.  Presumably this is 
  because the puromycin arm slightly disfavors annealing.

- The amount of non-specific annealing for o93 is pretty small and (I think) 
  unlikely to cause problems.  Given that I expect to get ≈90% efficiency for 
  the ligation reaction, 1% non-specific annealing should be negligible.

- I plan to use o93 for most of my ligation experiments, and to account for 
  non-specific annealing by running the appropriate controls and doing gel 
  densiometry.

