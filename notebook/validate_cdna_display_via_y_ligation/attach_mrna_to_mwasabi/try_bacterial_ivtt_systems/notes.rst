**************************
Try bacterial IVTT systems
**************************

2020/09/10:

So far, I've been using PURExpress for all of my cDNA display experiments.  
However, in light of the fact that PURExpress seems to be negatively affecting 
my mRNA (:expt:`66`), I think it's worth trying some other systems:

- PUREfrex: Another commercial PURE system sold by GeneFrontier.  According the 
  [Doerr2019]_ and [Niederholtmeyer2013]_, it has no detectable nuclease 
  activity (while PURExpress does have some).  It also has no His-tagged 
  components, so I can do His-tag purifications (but not reverse 
  purifications).

  Unfortunately, PUREfrex is not available in B2P, so ordering might be a 
  hassle.  I requested a quote.  See this page: 
  https://www.genefrontier.com/en/solutions/purefrex/lineup/

- Promega S30 extract, NEBExpress: Two commercial S30 lysates.  PURExpress 
  should have lower nuclease activity than any lysate, but maybe there's 
  something else going on.  Additionally, most mRNA display systems use lysates 
  (specifically, rabbit reticulocyte lysate), so it's not like lysates are 
  completely off the table.  Really I just have both on hand, so I feel like I 
  might as well try them.

I think it's worth trying rabbit reticulocyte lysate as well, but that would be 
a more involved experiment because I'd first have to reclone my constructs.

Results
=======

Observations:

- My f89 stock is not concentrated enough to reach the ideal concentration for 
  the NEBExpress reaction (400 nM).  The most I can reach is:

  .. math::

    \frac{\pu{1000 nM} \times \pu{1.44 µL}}{\pu{6 µL}} = \pu{240 nM}

  According to the results from :expt:`98`, that should be enough to get ≈50% 
  of the maximum expression.  Hopefully that will be good enough.

- I mistakenly entered the wrong target concentration for the PUREfrex 
  reaction.  As a result, I used 200 nM and not 500 nM.  I actually could not 
  have reached 500 nM, but I could've gotten to 350 nM.  Note that 200 nM 
  should still be enough to get 70% yield.
