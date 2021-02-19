***********************************
Check degradation via Northern blot
***********************************

One way to visualize mRNA is by Northern blotting.  Here I will try using this 
technique to determine whether the mRNA is full-length.

First try --- 2020/11/19
========================
.. protocol:: 20201119_optimize_crosslinking.txt

There wasn't any signal, but there are several possible explanations:

- The transfer may not have worked.

- The probe could've been quenched, because I didn't keep it in the dark.

Optimize transfer via TLC --- 2021/02/17
========================================
.. protocol:: 20210217_optimize_transfer.txt

.. datatable:: 20210217_optimize_transfer_notes.xlsx

Observations:

- The TLC method was capable of visualizing the bands (barely), but wasn't 
  sensitive enough to really judge if the transfer was complete.

Optimize transfer via Cy5 --- 2021/02/18
========================================
I realized that I could use Cy5 to very accurately monitor the transfer:

.. protocol:: 20210218_optimize_transfer.txt

  Because the gel tore, I wasn't able to directly compare the +/− transfer 
  bands to calculate what percentage of mRNA was transferred at each time 
  point.  To approximate this calculation, I took the following steps to 
  digitally "tear" the untorn band in the same way as the torn one.  This 
  should make the band intensities (roughly) comparable.

  - Open image in Gimp

  - Color > Auto > White Balance

    - To make the tear easier to see.

  - Select rectangle around tear, a bit wider than the band.

  - Switch to "select by color" tool, in "intersect" mode (Ctrl-Shift).

  - Select black color in rectangle; this is the tear.

  - Switch to "free select" tool, in "subtract" mode (Ctrl).

  - Remove any parts of the selection that aren't really part of the tear.

  - Make guide at right edge of selection

    - Note that the two bands aren't necessarily aligned, but this does help to 
      see how the tear should be positioned.

  - Copy the selection.

  - Undo the white-balance.

  - Paste the selection.

  - Drag selection so that it lines up about right with the untorn band.

  - Anchor.

  - Save as ``*_tear.tif``.

.. figure:: 20210218_optimize_transfer.svg

.. datatable:: 20210218_optimize_transfer.xlsx

  See caveats in protocol above; I had to manually account for a tear in the 
  gel.  These numbers should be considered qualitative.

Observations:

- The lower MW species transfer earlier.

- The transfer percentages should be taken with a grain of salt, because 
  different MW species transfer at different rates.  Because there is some 
  smearing, the transfer percentages don't exactly say how much of the *full 
  length* mRNA has transferred at a given time point.  

  Looking at the white-balanced images, it seems like all of the full-length 
  mRNA has transferred by 90 min, while some remains at 60 min.

- I tore the gel, which complicated analysis.  See the protocol section above 
  for details on how I accounted for the tear digitally when doing densiometry.

- My observations from yesterday (via the TLC method) are actually quite 
  consistent with these results.  I'm kinda surprised that I was able to get 
  roughly the right answer with such faint signal.

Conclusions:

- I think I'll do 60 min transfers going forward.  Even though I do think the 
  transfer is more complete at 90 min, I also think it's probably smart to keep 
  the transfer as short as possible.  The semidry transfer uses very little 
  buffer, so I assume that it's relatively prone to heating up (which I further 
  assume could cause problems in general).  This effect would not have been 
  captured by this experiment, because I basically took the gel out and cooled 
  it down every 30 min.  
  
  Most Western blotting semidry transfer protocols call for 15-30 min 
  transfers, but I was able to find some that call for 60 min.  Even though 
  Western and Northern transfers aren't really comparable (the buffers are very 
  different), this makes me skeptical of doing 90 min transfers.

  If I feel like I need the transfer to be more efficient at some point down 
  the line, I could either (i) try 90 min, (ii) do a wet transfer, or (iii) 
  load more material in the first place.  Note that wet transfers are generally 
  recommended Western blotting for large proteins.

