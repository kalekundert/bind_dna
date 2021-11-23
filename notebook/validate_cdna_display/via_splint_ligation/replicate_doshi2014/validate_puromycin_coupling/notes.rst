***************************
Validate puromycin coupling
***************************

My goal is to do the coupling reaction in the same way as [Doshi2014]_.  
Instead of purifying with the oligo-dT beads after purification, though, I used 
a desalting column to allow for −puromycin controls.

.. protocol:: 20211112_validate_puro.pdf 20211112_validate_puro.txt

.. figure:: 20211115_validate_puro.svg

- I don't see any evidence of coupling.  I also don't really see any evidence 
  of protein expression.

- The −linker control has a faint band around 33 kDa.  I don't really know what 
  this could be.

- The tRNA didn't run at the expected MW.  This is consistent with what I've 
  seen before (:expt:`117`), but I still don't know why.  The tRNA band is 
  sharper in tricine gels, so I might use that going forward.

- The Nb product and the unreacted tRNA run at almost the same MW.  This will 
  make it hard to visualize free Nb-GFP, although in principle it should still 
  be easy to see mRNA-conjugated Nb-GFP.


Next steps:

- Has my PURExpress gone bad?

  - Compare DHFR control for new/old PURExpress

- Can I successfully express Nb-GFP?

  - Titrate DNA and mRNA template.
  - Visualize with either:

    - FluoroTect + RNase
    - Western blot

- Does desalting cause any problems?

  - Optimized mRNA or DNA concentration.
  - +/− desalting column
  - Don't expect any problems, but I should do this control to be sure.
- 
