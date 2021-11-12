*****************
4 fragment Gibson
*****************

2021/10/29:

.. figure:: 20211029_check_f147_f153.svg

- f149 did not amplify well.  This is unexpected: f149 is basically the same 
  construct as f138, just with slightly different primer overhangs.  f138 
  amplified pretty cleanly with the same thermocycler conditions (see binder: 
  2021/09/21; see raw gel data: 2021/09/22).

  - I wonder if I just made a mistake setting up the reaction.  I'll set up a 
    thermocycler gradient and see what I get.  I'll also try setting up a 
    touchdown PCR at the same time.

  - It's not clear if the short band is just the primers, or a short 
    amplification product.

    According to `Thermo's Multiple Primer Analyzer`__, there is the 
    possibility of a primer dimer for this reaction::

      285_AMP_SR045_TM65 with 286_AMP_SR123_TM67
      285_AMP_SR045_TM65
      5-cggatcgaacttaggtagcccactcgataggtacaaccgg->
                                         | |||||
                                       <-gatggccactgcagccttaacttttcagggttactcacgg-5

    This same possibility exists with the o260/o261 primers used to amplify 
    f138, however.

    __ https://www.thermofisher.com/us/en/home/brands/thermo-scientific/molecular-biology/molecular-biology-learning-center/molecular-biology-resource-library/thermo-scientific-web-tools/multiple-primer-analyzer.html

- In contrast to f149, f147 and f148 amplified about the same as their f136 and 
  f137 counterparts.

- f150-f153 all amplified well.

2021/11/01:

.. protocol:: 20211101_pcr_gradient_touchdown.pdf 20211101_pcr_gradient_touchdown.txt

.. figure:: 20211101_check_f149.svg

- I repeated the amplification using gradient and touchdown PCR, and got good 
  amplification for every condition.  

  - One difference was that I setup the reactions on ice this time.  I'm going 
    to test to see if that makes a difference.

  - I could've also just messed something up last time.

2021/11/01:

.. protocol:: 20211101_pcr_ice_vs_rt.pdf 20211101_pcr_ice_vs_rt.txt

.. figure:: 20211101_f149_ice_vs_rt.svg

- Setting up the reaction on ice gives significantly more product than setting 
  it up at room temperature.  This is very consistent with primer dimers like 
  those predicted above being filled in by the polymerase before the reaction 
  starts.

- I didn't pre-heat the thermocycler before starting the reaction: that might 
  have given even cleaner product.

- In the future, I'm going to order hot-start Q5 and use it for everything.
