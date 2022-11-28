**********
Clone p239
**********

My p238 miraprep has significant RNA contamination.  That in itself doesn't 
really matter, since (i) RNA shouldn't interfere with cloning reactions and 
(ii) I need to do a gel purification after the restriction digest anyways.  But 
this does make it hard to know how much restriction enzyme to use.  

I decided to measure the DNA concentration by running a gel.  This gave me a 
concentration of 15 ng/µL, which is kind of hard to believe.  Is it possible 
that the fast-running RNA contamination depleted the gel of dye, such that the 
later bands are artificially dim?  No, that would make the dilutions 
non-linear, which they aren't.  By the same measurement, the fast-running band 
is only ≈50 ng/µL (data not shown), although if it's RNA the it wouldn't be 
expected to bind the dye as well.  Still, that's nowhere near the 5000 ng/µL 
that the nanodrop measured.

Something's wrong.  I don't fully trust the nanodrop or the gel.

I'll setup a 50 µL digestion reaction.  By the gel, that would be 0.75 µg.  
I'll maybe use enough enzyme to digest 10 µg (250 ng/µL at 40 µL), though, to 
be safe.

2022/09/23
==========
- I've purified both f190 and f191, but both have very low concentrations.  As 
  a result, I'll be using about 10x less DNA than recommended for this 
  transformation.  That should still be enough for this small library, but it's 
  something I'll want to improve on.

- Regarding the insert (f191):

  - I can always do more cycles of PCR if I need more product.

  - The digestions reaction didn't seem to go to completion.  I'm experimenting 
    with that now: :expt:`198`.

  - I didn't measure yield for any of the purification steps.  Doing that would 
    be the first troubleshooting step.

- Regarding the backbone (f190):

  - The miniprep just doesn't give very good yield, see :expt:`194`.  The only 
    remedy I can think of is to switch to a different plasmid backbone.

2022/09/27
==========
.. protocol:: 20220926_make_p239.pdf 20220926_make_p239.txt

.. datatable:: 20220927_electrotransform_p239.csv

- Purifying the insert may improve transformation efficiency, if I can recover 
  more DNA.

  - For this transformation, I PAGE-purified the insert.

  - The yield from the purification was poor, so I transformed ≈10x less DNA 
    than the minimum recommendation.

  - Despite this, I still got a comparable number of colonies to the p237 
    transformation (for which I used crude insert).  This suggests that the 
    purified insert gave rise to a more efficient transformation.

    - The cells I used for this transformation had already been thawed once, so 
      that could be the cause for a slight decrease in efficiency.

    - It might be worth comparing the quality of the two libraries.

    - On the other hand, one could argue that transformation efficiency doesn't 
      really matter for these libraries, so I might as well skip the extra 
      purification steps and just use crude insert.

- My midiprep yield was poor, but probably good enough.

  - Yield: 1.9 µg; 20 µL at 95 ng/µL
  - The 260/280 ratio was good: 1.83
  - I used the gravity column, but I think the plunger columns might be easier.
  - Qiagen recommends taking aliquots throughout the process to troubleshoot 
    issues like this.  I didn't take any aliquots this time, but it's porbably 
    worth doing next time.

2022/10/19
==========
.. protocol:: 20221019_make_f190.txt

.. figure:: 20221019_gel_purify_f186_f190.svg

- I had two tubes of p238 in my freezer.  The one I used was ≈200 ng/µL with a 
  poor 260/280 ratio, and apparently had very actual plasmid (see above gel).  
  The other was a miraprep.  In the future, I should make sure I have a good 
  plasmid prep before doing this digest.

- It seems like the digestion went to completion, although I (inadvertently) 
  used a massive excess of enzyme.

.. protocol:: 20221020_make_p239.txt

.. datatable:: 20221021_electrotransform_p239.csv

- I manually counted 1701 colonies, but calculated (from my dilutions) that 
  there should be 3125 colonies.  I'm not sure why there's a discrepancy...

  - I checked by hand that my calculations are in the right ballpark, so it's 
    not that there's a bug in my script.

- I got 10x more transformants with insert than without.

  - The −insert count is probably not very accurate, though, because I only got 
    3 colonies total.

  - I'm not sure if this is good or not.  Clearly it's better than getting 
    ≈equal colonies, but 10% of the library missing the insert still seems like 
    a lot.

- Even though I didn't get many colonies, the transformation efficiency wasn't 
  horrible:

  - I used 3.23 µL × 0.288 ng/µL = 0.93 ng backbone.
  - I got 1.7×10⁴ transformants (accounting for the fact that I only plated 100 
    µL).
  - That works out to 2×10⁷ transformants/µg DNA.
  - NEB advertises 2×10¹⁰ transformants/µg DNA, when using 1 ng pUC19.
  - I'm 1000x worse that that, but I'm using ligated plasmid, so perhaps that's 
    around what would be expected.

- I'm confident I'll get better results if I just manage to purify more DNA.

2022/10/31
==========
Do a test digest to try to understand why my sequencing reactions all failed.

.. protocol:: 20221031_test_digest.txt

.. figure:: 20221031_test_digest_p239.svg

- There's probably at least some p239 in these reactions.

  - All of the reactions have at least a faint band at the expected MW, and 
    this band is not present in the negative control.

  - For HaeII, this band is quite strong.

- I think I might have the pUC backbone, somehow.

  - The only restriction enzyme that appears to cut all of the starting 
    material is HaeII.  This is also the only site that is present in p237.  
    It's present in the pUC ORI, so any plasmid with that backbone would be cut 
    by HaeII.

  - HaeII gives the expected 3.4 kb band, but also an unexpected 2.4 kb band 
    (and a number of very faint smaller bands).  The 2.4 kb band is consistent 
    (in length) with p237, and probably other plasmids of mine with the pUC 
    backbone.

  - On 10/19, I purified f186 (pUC backbone) and f190 (pSC101 backbone) at the 
    same time.  Perhaps they somehow contaminated each other, and I should try 
    purifying f190 again...

- I don't think there's any problem with p238 (which is where the backbone 
  comes from).

  - I have Sanger sequencing data but not full-plasmid Nanopore data for p238.  
    The Sanger reads overlap the Rep101 ORF by ≈100 bp, which makes it hard to 
    believe that p238 has the wrong plasmid backbone.

2022/11/08
==========
.. protocol:: 20221108_make_f190.pdf 20221108_make_f190.txt
.. figure:: 20221108_check_digestion_f190.svg
.. figure:: 20221108_gel_purify_f190.svg

- The plasmid seems to be completely digested.

  - The original circular bands are completely gone, indicating that the 
    plasmid is completely linearized.

  - The 0.4 kb band is clearly visible, and there don't seem to be both 3.6 and 
    3.2 kb bands.  This suggests that the dominant species has been cut twice.

- The gel purification has a bit of a smear, but the band should be 
  predominantly the expected species.

2022/11/09
==========
.. protocol:: 20221109_make_p239.pdf 20221109_make_p239.txt
.. datatable:: 20221110_electrotransform_p239.csv

- Very few transformants, presumably because the insert was so dilute.

2022/11/14
==========
The Sanger sequencing failed, but I decided to continue trying to make f200 and 
then p243 anyways.  See :expt:`204`.  

I bet my sequencing failed because the template concentration was too low.  The 
smearing visible in the gel suggests that the concentration measured by the 
Nanodrop was an overestimate.  Maybe in the future I'll want to try PCR 

2022/11/15
==========
I transformed my leftover ligation reactions, in the hopes that a midiprep 
would give better DNA somehow:

.. protocol:: 20221114_retransform_p239.pdf

.. datatable:: 20221115_electrotransform_p239.csv

- Lyophilizing the DNA before doing an electrotransformation seems to be a bad 
  idea.

  - The decay times reported by the electroporator were short, and the 
    transformation efficiency was worse than if I'd just used 2 µL without 
    concentrating.

  - Presumably this is because residual salts in the sample were concentrated 
    enough to cause problems, although I eluted in water, so I'm surprised that 
    there could be enough residual salts to matter.

- The sequencing data for this prep is messy, but seems to be the right thing.  
  I think the problem I've been having is just that this is a low copy plasmid, 
  and I need enough to sequence.

- Next time, I should do the "low-copy" midiprep protocol, which I think just 
  uses more cells and more buffers.
