**********************
Measure standard curve
**********************

My ultimate goal is to develop a high-throughput assay that will give me 
accurate measurement of protein binding to DNA.  In order to know how accurate 
such an assay is, I'll need to measure standard curves.  My goal for this 
experiment is to measure standard curves for all of the different assay 
variants I have in mind (e.g. different reporters, different plasmid 
refactorings, etc.) and compare which gives the largest and most linear dynamic 
range.

Considerations
==============

Published data
--------------
I'd rather build my standard curve from existing binding data, rather than 
having to collect my own.  These references vary both the target site and the 
protein itself.  I'd probably prefer to vary just the target site, because that 
would be easier and cheaper to clone, but in the long run I'll probably want to 
do both:

- [Jamieson1994]_

  .. datatable:: jamieson1994_zif268.xlsx

  - :math:`K_D` measured by gel shift.
  - These data were used by [Jamieson1996]_ for calibration, which is basically 
    what I want to do as well.

- [Wu1995]_, [Yang1995]_

  .. datatable:: wu1995_zif268.xlsx

  - :math:`K_D` measured by SPR.  To my knowledge, this is the most accurate 
    way to measure :math:`K_D` [Matos2010]_ [Yang2016]_.
  - Table 2 in [Wu1995]_ includes all the data from [Yang1995]_
  - Wildtype measured with 4 different targets, but all 3 non-WT targets have 
    similar :math:`K_D`.

- [ElrodErickson1999]_

  .. datatable:: elroderickson1999_zif268.xlsx

  - :math:`K_D` measured by gel shift.
  - Measuring :math:`K_D` is the focus of the paper.
  - Alanine mutant of base-contacting residues.
  - Only 3 target sites, but they span a good range.
  - Has references to 7 other papers that measured :math:`K_D` for Zif268, and 
    claims that the measurements for wiltype protein/target range from 0.01-6.5 
    nM (depending on the specific conditions).

- [Bulyk2001]_

  .. datatable:: bulyk2001_zif268.xlsx

  - :math:`K_D` measured by phage ELISA.
  - I probably trust these data less than the gel shift data from 
    [Jamieson1994]_, since in this experiment the Zif268 is still tethered to 
    phage.
  - But these data also span a much larger range.

Which target sites to use?

- Kinda want a geometrically increasing series of binding constants, e.g. 1 nM, 
  2 nM, 4 nM, 8 nM, 16 nM, etc.

- ``gcg tgg GAG``:
    
  - Measured by [ElrodErickson1999]_ (gel shift) and [Jamieson1994]_ (gel 
    shift).  Somewhat different :math:`K_D` relative to wildtype (3x vs 10x), 
    but in both cases this is very nearly the second best target.  So this 
    sequence is probably a good test of an assay's ability to distinguish 
    between two pretty similar binders.

- ``gcg tgg GCC``:
  
  - Measured by [ElrodErickson1999]_ and [Jamieson1994]_, but they disagree 
    significantly about how well it binds.

- ``gcg TTG gcg``:

  - Measured by [Wu1995]_ (SPR, 83.7 nM, 13x worse than WT) and [Bulyk2001]_ 
    (EMSA, 71 nM, 23x worse than WT).  Here the agreement is quite good.  That 
    makes me more inclined to trust the other [Bulyk2001]_ numbers, even though 
    they're measured while still bound to phage.

- No other targets were measured by multiple references (that I found).

- Based on what I have so far, I'm inclined to use:

  - ``gcg TGG gcg``: :math:`K_D` = 3.0 nM
  - ``gcg TAG gcg``: :math:`K_D` = 6.7 nM
  - ``gcg TTG gcg``: :math:`K_D` = 71 nM
  - ``gcg CAG gcg``: :math:`K_D` = 380 nM

  - These are all from the [Bulyk2001]_ EMSA measurements.  ``TTG`` is also 
    measured in [Wu1995]_.
  - I'm more inclined to trust [Bulyk2001]_ because the one measurement that 
    overlaps with the SPR data has good agreement.
  - 

- Which protein variants to use?

  - [Wu1995]_, which I think is the highest quality :math:`K_D` data, has a 
    handful of variants with a good range of activities:

    - WT
    - C7
    - C9
    - C10
    - F8
    - F15

    - I could use these variants with both the wildtype ``gcg tgg gcg`` target 
      site and the ``gcg tgg TGT`` site.


