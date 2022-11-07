***********************
Compare competent cells
***********************

2022/04/19:

I made a batch of Mix & Go competent MACH1 cells today, and I wanted to compare 
their competency to the two other MACH1 competent cell preps I have in my 
−80°C:

- MACH1 cells from Thermo that I received in 2019.
- CaCl₂ competent MACH1 cells that I made on 2022/02/28.

Part of my goal for this experiment is to establish a baseline, so in the 
future I can repeat the experiment and see how long the cells remain competent.

.. protocol::

  - Prewarm 3 LB+Carb plates at 37°C for >1h.

  - Thaw the following volumes of competent cells:

    - Mix & Go: 25 µL (based on Tina's recommendation)
    - Thermo: 10 µL (based on what I usually do)
    - CaCl₂: 100 µL

  - Add 1 µL p2 (230 ng/µL pUC19) to each transformation.

  - Incubate on ice for 5 min.

  - Heat shock (42°C for 45 sec) the Thermo and CaCl₂ transformations.

  - Add 100 µL SOC to the Thermo transformation.

  - Plate the transformations.

    - Note that I normally add SOC to the CaCl₂ transformation, but didn't this 
      time.

.. figure:: 20220420_compare_comp_cells_mach1.tif
.. figure:: 20220420_compare_comp_cells_cacl2.tif
.. figure:: 20220420_compare_comp_cells_mix_n_go.tif

Observations:

- The Thermo cells are barely competent anymore, but that's expected because 
  they're so old.  It's hard to see in the above image, but I actually got 2 
  colonies for that transformation.

- Right after I made the CaCl₂ cells, they were more competent than the Thermo 
  cells (I didn't save any data about this, but I transformed the same cloning 
  reaction with both and got more colonies with the former).  Now the CaCl₂ 
  cells seem almost completely inactive.  I didn't get any colonies for this 
  transformation.

- The Mix & Go cells are clearly more competent than the alternatives.  I 
  counted 421 colonies, which is: :math:`\frac{\pu{421 cfu}}{\pu{230 ng}} 
  \times \frac{\pu{1000 ng}}{\pu{1 µg}} = \pu{1830 cfu/µg}`.  That's far less 
  than the 10⁹ advertised by Zymo, though.

