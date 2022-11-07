*********************
Compare miniprep kits
*********************

My minipreps recently have all had very large amount of RNA concentration, 
judging by:

- The 260/280 ratio measured by the Nanodrop.
- The actual concentration of plasmid, measured by gel densiometry.

To try to troubleshoot this, I want to miniprep the same culture with several 
different kits.  I'm going to borrow kits from Tina for this experiment, as her 
minipreps have been working well recently.

I also want to try miniprepping a smaller volume of cells, and the Qiagen 
troubleshooting webpage suggested that a smaller volume of cells would have 
less RNA, and would therefore be less likely to overwhelm the RNase in the 
miniprep protocol.

Results
=======

2022/09/14
----------
.. protocol:: 20220914_compare_minipreps.pdf

.. figure:: 20220914_compare_miniprep_kits_40ms.svg

.. datatable:: 20220914_compare_miniprep_kits_40ms.xlsx
   :sheet: Report

- The gel signal seems to get saturated at about 9000 px.  The 5 µL for a 
  number of the minipreps hits this threshold, while the 0.5 µL load only hits 
  it for 3131.

  - I took two images of the gel, with exposure times of 122 ms and 40 ms, 
    respectively.

  - The former has saturated pixels; the latter doesn't.

  - I double-checked that I used to 40 ms image (no saturated pixels) for the 
    densiometry analysis.

  - Despite this, the signal seems to saturate at peak areas of ≈9000 px.

  - The 5 µL and 0.5 µL loads pretty much agree for the most dilute preps: the 
    1 mL Qiagen preps and the Thermo prep.  These bands all have areas less 
    than 8000px.

    - To be honest, even the 7500 px bands may be a bit saturated.  They give 
      concentrations of 30 ng/µL, compared with 36 ng/µL for the 10x dilution.

  - I'm not sure what the reason for this is, but it goes to show that I should 
    always run a 10x dilution if I want to use a gel to quantify DNA 
    concentration.

- My standard curve is ok.

  - I ran the ladder in 3 lanes with a 16x range in concentration in the hopes 
    of getting a broad standard curve.  That didn't really work because:

    - The gel was too short, so the ladder bands weren't able to separate very 
      much.

    - The 2 µg lane ran weirdly, and wasn't much brighter than the 500 ng lane.

  - The standards range from 200-1300 px.  In contrast, my unknowns (that don't 
    appear to be saturated) are mostly in the range 900-5000.  Ideally, there 
    would be better overlap between these ranges.  I think it will be hard to 
    achieve that with ladders having so many bands, though.  Ideally I'd have a 
    1-3 bands of known concentration, but I haven't been able to find a product 
    like that (and I don't trust myself to make it).

  - All that said, the curve is pretty linear in the range it covers, and goes 
    pretty close to the origin.

- The two plasmids that Tina miniprepped (3131 and 1124) don't appear to have 
  any RNA contamination:

  - The gel densiometry concentrations agree with the Nanodrop concentrations.  
    3131 actually seems to be saturated in the gel, but the result is still 
    consistent with the Nanodrop measurement.

  - The 260/280 ratios are below 2.0.

- The problem isn't with my Qiagen kit.

  - I got pretty much identical results with Tina's Qiagen kit.

- The problem isn't with the number of cells I'm miniprepping.

  - Using 1 mL of cells instead of 5 mL:

    - I get 1/4 the yield.
    - ≈33% of the Nanodrop signal is accounted for by the gel, instead of ≈25%.  
      This is consistent with a slight increase in purity.

  - The decrease is yield is nearly proportional to the decrease in cells.
  - The slight increase in purity is consistent with the RNase being saturated, 
    but even with just 1 mL cells, the purity is not good.

- The 100 ng/µL yield from my Qiagen prep isn't *that* bad, although it's still 
  less than I would expect, and it does seem the be 75% RNA.

  - This is also 2-10x better than I've seen previously.

  - Given that my 5 µL load seem to be saturating something, and this is the 
    first time I included dilutions of the unknowns, it's possible that some of 
    my previous gel densiometry measurements have been inaccurate.

  - That said, in all my previous experiments, I only loaded 1 µL of DNA, which 
    wouldn't have really been saturating in this experiment.

- The NEB prep is pretty clean, but still low-yield.

  - The low yield may be because I only used 1.5 mL of cells.  That is what the 
    protocol recommended, and 5 mL would have been over the explicit maximum in 
    the protocol, but I think Tina uses 5 mL with this kit.

  - The 260/280 ratio suggests that the prep is not clean, but I put more trust 
    in the fact that the Nanodrop and the gel agree pretty well on the 
    concentration.

- I don't understand how Tina can routinely get 900 ng/µL from a miniprep.

  - `Qiagen miniprep columns`__ have a capacity of <20 µg.  Assuming a 50 µL 
    elution volume, that's 400 ng/µL.

    __ https://www.qiagen.com/us/products/discovery-and-translational-research/lab-essentials/plastics/qiaprep-spin-miniprep-columns/

- My plasmids do seem to have contamination.

2022/09/23
----------
.. protocol:: 20220923_compare_miniprep_qiagen_zymo.pdf

.. figure:: 20220923_compare_minipreps_fits.svg

.. datatable:: 20220923_compare_miniprep_qiagen_zymo_report.xlsx

- I don't trust the densiometry data for this experiment:

  - The standard curve has a much lower R² coefficient than usual.

    - It's notable that only the top-most bands (3-10 kb) in the first ladder 
      lane are saturated.  All of the ladder bands have similar amounts of DNA 
      (≈40 ng for the normal bands, ≈120 ng for the bright ones), so if the 
      top-most bands are saturated, the bottom-most bands should be too.

  - The dilutions of the Zymo prep aren't self-consistent.

  - Maybe I didn't soak the gel for long enough?

- The Zymo kit gives less yield than the Qiagen kit (by Qubit), but seems much 
  cleaner (by Nanodrop).

  - The contamination (whatever it is) must be due to the Qiagen kit, since 
    using a different kit gives clean product.

  - But the low yield is most likely intrinsic to the plasmid.  That's 
    surprising to me, since this is frickin' pUC19, but I really can't think of 
    any other explanation at this point.
