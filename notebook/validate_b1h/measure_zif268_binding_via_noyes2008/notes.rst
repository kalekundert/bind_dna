***************************************
Measure Zif268 binding via [Noyes2008]_
***************************************

My goal here is to simply recapitulate the basic B1H assay using the controls 
from [Noyes2008]_, to make sure the strain and the plasmids work in my hands.  
To that end, I'm planning to individually transform and test target and 
non-target bait plasmids.

2021/11/16:

.. protocol:: 20211019_make_nm_agar.txt 20211102_make_nm_media.txt 20211104_test_controls.pdf 20211104_test_controls.txt

.. figure:: 20211109_noyes2008_controls.svg

  s4: USE hisB- pyrF- rpoZ- target=TGG; s5: USE hisB- pyrF- rpoZ- target=AAA.  
  I counted colonies from the images, not from the actual plates.  This is 
  probably less accurate, and in some cases there were colonies obscured by the 
  grid lines that I couldn't count.

Observations:

- There were random colonies growing all over the plate, so I should take small 
  counts with a grain of salt.  I'm not sure what I could've done to prevent 
  this.

- The TGG-target strain grew just as well in the 4 mM 3-AT condition as in the 
  +His condition.

  I think this is expected.  Looking at the [Noyes2008]_ data (Fig S5), the 
  positive control is pretty much unaffected by 3-AT, even up to 50 mM (0 mM: 
  267 colonies; 50 mM: 235 colonies).  
  
  I don't remember why I chose to use 1-4 mM 3-AT concentrations, but I think 
  it was to measure the sensitivity of the AAA-target strain.  For reference, 
  the GCG-, GAG-, and GCA-target strains (3rd finger) are all mostly 
  insensitive to 50 mM 3-AT, while the GGG-target strain is sensitive to 25 mM 
  3-AT ([Noyes2008]_, Fig S5).  The AAA-target strain should be much weaker 
  than any of these variants, hence the much lower 3-AT range.

- The AAA-target strain did not grow in the +His condition.

  I think this is because there was no uracil in the plates.  The NM −His −Ura 
  plates are selective for both HIS3 and URA3, so my "permissive" plate was 
  still selective by virtue of being −Ura.  

  The question is how to deal with this.  [Noyes2008]_ includes 200 µM uracil 
  in the liquid resuspensions before plating the cells, which means that some 
  poorly-defined quantity of uracil would be present in all of the plates.  
  (The protocol is also unclear about diluent used for the final dilution of 
  the cells.)  Overall, I don't think the [Noyes2008]_ protocol is very good in 
  this regard.

  A reasonable alternative is to include uracil in the "permissive" plate.  So 
  the "selective" plates have different stringencies based on [3-AT], and the 
  permissive plate is really permissive.  I already know that the AAA-target 
  variant won't grow in the absence of uracil, even with 0 mM 3-AT, so while I 
  should repeat the experiment, I already know what will happen.  It might be 
  good to include the GGG-target in this experiment, although I'd if I do that 
  I'll need to include a wider range of 3-AT concentrations to see a titration.

  I don't like the idea of including any uracil in the selective plates, 
  because the library selection acts on both His and Ura.  I think this 
  experiment should mimic that.

- I had to incubate the plates for 36h to see colonies.  It turns out that 
  [Noyes2008]_ also calls for a 36-48h incubation; I just didn't copy that 
  detail into my protocol.

Takeaways:

- If I need less selective conditions, I could try +Ura −His +3-AT.  I don't 
  see an immediate need for this, though.

Next steps:

- Repeat experiment with uracil in permissive plates.

- Measure OD.

  [Noyes2008]_ measures OD in addition to counting colonies, to quantify colony 
  size.  This seems like a reasonable idea, but it requires that condition gets 
  it's own plate.  I think that for my purposes, which are just to show that 
  the assay works qualitatively, I won't get much benefit from doing this.

- Make exact same media as [Noyes2008]_.

  - I cut a corner by use a −His −Ura amino acid mix.
  - The amino acid mix used by [Noyes2008]_ includes more amino acids.
  - I don't think this is really the problem, though.

- Clone a −Zif268 negative control (i.e. untethered ω-subunit).  [Noyes2008]_ 
  used this as a negative control in their library selections.

- Clone the Zif268 targets used in [Noyes2008]_, i.e. GGG.

- Include single-plasmid constructs in the experiment.

  - What can go wrong?

