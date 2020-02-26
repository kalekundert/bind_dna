***************************
Express mRNA via PURExpress
***************************

In :expt:`express_zif268_via_naimudden2016` (2020/02/18), I saw poor expression 
of protein from mRNA.  Here I want to determine how much mRNA is needed to get 
appreciable expression with PURExpress.  I decided to start with mWasabi, 
because it should be pretty easy to detect its expression.  I even think it 
might remain fluorescent in after SDS-PAGE based on this `high school/undergrad 
lab kit 
<https://www.bio-rad.com/en-us/product/pglo-sds-page-extension?ID=a41608e9-b348-43e0-98bb-d0ae12664e06>`__, 
which would be extra convenient.  Even if that doesn't work, though, it's still 
a protein that I know to be well-expressed.  DHFR (the NEB control) is another 
option.
  
Considerations
==============
NEB recommends using 1--5 µg mRNA per 25 µL PURExpress reaction.  This 
corresponds to 1.6 µL of 10 µM f11 per 10 µL reaction.

Results
=======

2020/02/25
----------
.. protocol:: 20200225_split_purex_page.txt

.. figure:: 20200225_express_mwasabi_mrna_10_10.svg

   PURExpress: Steps 2-3 in the above protocol.

- mWasabi remains fluorescent in SDS-PAGE gels.

- mWasabi runs a little higher than it should, but I'm not going to read 
  anything into that.

- It's strange that the 1 µM reaction appears to have produced significantly 
  more mWasabi than the 10 µM reaction.  Some explanations:

   - I mixed up the reactions (almost certain that this didn't happen)

   - SDS is somehow interfering with mWasabi fluorescence such that the 
     fluorescent signal is not actually proportional to the amount of protein 
     present.  The eliminate this possibility, I'm staining the gel with 
     Coomassie.

   - Using less mRNA gave the proteins more time to fold correctly, so I ended 
     up with more functional protein.  I'd want to repeat this experiment to 
     verify this, with more of a concentration gradient.

     .. todo::

         Repeat the mWasabi mRNA expression experiment with more mRNA 
         concentrations.

  NEB does say that the optimal amount of mRNA needs to be optimized for each 
  gene, so maybe that's just what I'm learning the hard way here.

2020/02/26
----------
I'm going to do a serial dilution of mRNA from 10 µM to 0.1 µM in 7 steps.  I 
know that expression is higher at 1 µM than at 10 µM, so I'm hoping that 0.1 µM 
is low enough to see expression go back down again (since I want to find the 
maximum).   I chose 7 steps for several reasons:

- To get reasonably fine grained data.

- To include 1 µM, to compare to yesterday's experiment.

- To use exactly 2 aliquots of PURExpress (with 5 µL reactions).

.. protocol:: 20200225_serial_purex_page.txt
