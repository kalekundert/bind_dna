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
