******************************
Check degradation via broccoli
******************************

One way to visualize mRNA is the incorporate a broccoli aptamer and the stain 
with DFHBI-1T [Filonov2015]_.  Here I will try using this approach to determine 
whether or not the mRNA is full length.

Considerations
==============

Broccoli vs. tBroccoli
----------------------
[Filonov2015]_ uses both broccoli and tBroccoli (broccoli fused to a tRNA 
scaffold to help improve folding).  I decided to use plain broccoli because 
[Filonov2015]_ shows that the tRNA scaffold is recognized by certain RNases in 
bacteria (specifically RNases E, T, and PH).  The broccoli aptamer seems to 
fold just as robustly, and it's shorter.

Sequence: AGACGGTCGGGTCCAGATATTCGTATCTGTCGAGTAGAGTGTGGGCT

Results
=======

Visualize mRNA via broccoli --- 2020/11/09
------------------------------------------

.. protocol:: 20201110_check_degradation_via_broccoli.txt

.. figure:: 20201109_confirm_dfhbi_stain.svg

.. datatable:: 20201109_confirm_dfhbi_stain.xlsx

- The DFHBI-1T staining worked.

- I got much lower sensitivity than [Filonov2015]_.

  The smallest quantity of mRNA I could detect on this gel was 8 ng.  
  [Filonov2015]_ claims that quantities as low as 100 pg can be detected.  This 
  implies that my signal is about 80x less than it should be.

  One likely problem is that I used Tris instead of HEPES for the staining 
  buffer.  (I just ran out of time to make the HEPES buffer.)  Tris is known to 
  chelate magnesium and many other metal ions [Fischer1979]_.  A quote from the 
  abstract: "great reservations should be exercised in employing Tris as a 
  buffer in systems which also contain metal ions".  This may have had the 
  effect of reducing the amount of magnesium available to the broccoli aptamer, 
  preventing it from folding correctly.

Check degradation --- 2020/11/10
--------------------------------

.. protocol:: 20201110_check_degradation_via_broccoli.txt

.. figure:: 20201110_check_degradation_via_broccoli.svg

Observations:

- The signal in the GFP/DFHBI-1T channel is very faint and greatly amplified in 
  the image shown.

- The signal in the Cy5 channel is also amplified to make the faint bands in 
  the +mRNA lanes visible.

- The band representing the mRNA-linker conjugate (≈198 kDa) is very faint.

Conclusions:

- I think the broccoli tag is interfering with protein expression.  I don't see 
  any evidence of GFP expression, which I've seen easily in similar 
  experiments, e.g. :expt:`65`.

  - I entered my mWasabi mRNA sequences (with and without broccoli) into 
    RBScalc, and both were predicted to be efficiently translated.  In fact, 
    the broccoli+mWasabi sequence was predicted to be better expressed.

  - It might be worth repeating this experiment with f85 as a positive control.

- The linker ligation didn't work well.  In many previous experiments, I've 
  seen nearly all of the linker be attached to the mRNA.  In this case, only a 
  tiny fraction was attached.

  - Maybe I lost some material using the spin column.  I can probably skip that 
    step for now, since there's very little unreacted linker in my reactions to 
    begin with.

    - This can't explain why there was no GFP expression in the -linker lanes, 
      though, because I'm 95% sure I didn't use the spin column on that mRNA (I 
      did this experiment a week ago, so I can't perfectly remember what I 
      did.)

  - It's strange that the mRNA appears in the DFHBI-1T channel in the −linker 
    lanes, but not the +linker lanes.  Based on the faintness of the Cy5 
    signal, I'm guessing that the quantity of mRNA in the +linker lanes is 
    simply below the detection limit of the DFHBI-1T stain.

Check degradation --- 2020/11/16
--------------------------------
Repeat the above experiment with mWasabi mRNA (−broccoli) as a positive 
control:

.. protocol:: 20201116_check_degradation_via_broccoli_v2.txt

.. figure:: 20201116_check_degradation_via_broccoli.svg

Observations:

- Both channels are amplified.

- The mRNA without the broccoli aptamer behaves as expected:

  - The poly-A linker is ligated successfully.
  - mWasabi is expressed successfully.

- The mRNA with the broccoli aptamer behaves as it did previously:

  - No apparent ligation with the poly-A linker.
  - No apparent expression of mWasabi.

- The linker appears less homogeneous in PURExpress than in water.

- Some of the DNA ladder bands are hard to see in this image, but they all show 
  up clearly when the blue channel is made brighter.

Conclusions:

- The −broccoli positive control confirms that nothing was wrong with the 
  ligation or expression reactions.  This suggests that the broccoli aptamer is 
  somehow inhibiting both ligation and expression:

  - It's easy to rationalize how the broccoli aptamer could interfere with 
    expression, since it's adjacent to the RBS.

  - I can't think of any reasons why the aptamer would interfere with ligation, 
    though.  The ligation site is on the completely opposite side of the mRNA.

- It's interesting to compare the Cy5 bands in the −PURExpress lanes with and 
  without broccoli.  The unligated band (≈15 kDa) is about the same intensity 
  in both lanes, but the ligated band (≈10 kb) is much brighter in the 
  −broccoli lane.

  - There must be about the same amount of mRNA and linker in both lanes, 
    because there were no purification steps.

  - I'm tempted to say that something to do with FRET or quenching is going on, 
    but I've used GFP with Cy5 before without having this problem, e.g.  
    :expt:`36`.  Maybe the fluorophores were in closer proximity in this 
    experiment, though.

  - This would be easy to test: just image the Cy5 channel before staining.

- I didn't image the gel before staining with DFHBI-1T because I was running 
  short on time, but I really wish I had.  That would have told me if the lack 
  of Cy5 in the +broccoli, +linker, −PURExpress lane is due to some sort of 
  quenching, or the linker really not being there.

- I should've used to ssRNA ladder instead of the dsDNA ladder...

Discussion
==========
I'm going to give up on using broccoli to visualize the mRNA.  I don't know 
exactly what the problem is, but I clearly have reason to doubt if results from 
+broccoli mRNAs could be applied to −broccoli mRNAs.
