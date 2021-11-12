******************************
Remove unligated mRNA via PAGE
******************************

2021/05/04:

Given the trouble I've been having with the puromycin coupling reaction, I 
think it might be helpful to start with mRNA that is 100% ligated to the 
linker, e.g. no free mRNA and no free linker.  

In this experiment I'm going to look at PAGE purification, which has been used 
by [Reyes2021]_ and others for this exact purpose.

Measure yield
=============
2021/05/19:

I wanted to begin by gel extracting pure mRNA, to practice the protocol and to 
get a sense for the yield I can expect.

.. protocol:: 20210519_measure_yield.pdf 20210519_measure_yield.txt

.. datatable:: 20210519_measure_yield.xlsx

Observations:

- 60-75% of the input mRNA was recovered.

  That's much better than I was expecting.  For the ligation reaction, I can 
  probably expect a total yield of 40% (ligation reaction) × 70% (gel 
  purification) ≈ 30%.  That should be workable.

- The smaller mRNA (FLAG) was recovered slightly more efficiently than the 
  larger mRNA (Zif268).  Some thoughts:

  - This may just be because the smaller RNA diffuses out of the gel more 
    easily.

  - This may also be because the f111 gel slice was smaller (the band was less 
    diffuse, presumably because less mass was loaded).  I didn't weigh the gel 
    pieces, but maybe equilibrium favored more mRNA leaving the gel because the 
    solvent:gel ratio was higher.

- The gel was not overrun with 10 µM mRNA, which is about as concentrated as I 
  can imagine using.

Gel purify f112 and f118
========================
.. protocol:: 20210526_make.pdf 20210526_make.txt

.. protocol:: 20210527_gel_purify_f112_f118.pdf 20210527_gel_purify_f112_f118.txt

.. figure:: 20210610_gel_purify_f112_f118.svg

.. datatable:: 20210527_gel_purify_f112_f118_yield.xlsx

Observations:

- The purity is pretty good for the Zif268 mRNAs.  I was worried about these, 
  because the bands didn't have a lot of separation, but 90% is probably good 
  enough.  I could also try running the gel longer to get better separation.

- The yield was slightly lower than I expected.  I expected to get ≈40% yield 
  for the ligation step and ≈70% yield for the PAGE purification step, for an 
  overall yield of ≈30%.  I got close to that for the FLAG constructs (20-25%), 
  but not as close for the Zif268 constructs (10%).

  I don't know whether the ligation of the PAGE purification is more 
  responsible for the low yield.  To know in the future, I'd have to run a 
  lane with a small amount of each ligation reaction.  That would allow me to 
  calculate exactly how much mRNA was ligated, and from that I could back out 
  how efficient the PAGE purification was.

- I ended up loading:

  - 25 µg f11
  - 8 µg f111

  This was a little more than I planned to load (I'd read somewhere to not 
  exceed 20 µg/lane).  The gel didn't seem overloaded or smeared, though.   
  I'll have to run another gel on the purified products to see how well the 
  purification worked.  If I scale up the purification, I'll probably want to 
  split the samples between multiple lanes.

Conclusions:

- PAGE purification seems like a viable way to remove unreacted mRNA.  I'll 
  probably want to scale up the reaction a little bit more (e.g. 2x) next time, 
  but I'm able to get reasonably concentrated mRNA.
