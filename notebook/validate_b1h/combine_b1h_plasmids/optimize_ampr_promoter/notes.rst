************************
Optimize AmpR expression
************************

In :expt:`158`, I suspect that some of the results can be explained by the AmpR 
promoter being too weak in the pSC101 backbone.  Here I want to consider this 
possibility in more detail and design constructs that address this problem.

Copy number
===========
My AmpR gene comes from a pUC vector, which has a copy number of ≈500 according 
to Addgene_.  pSC101 has a copy number of ≈5 according to the same source, 
suggesting that I might aim to increase expression by 100x.

.. _Addgene: https://blog.addgene.org/plasmid-101-origin-of-replication

Existing/buffered promoters
===========================
I copied this AmpR gene from a high-copy vector (pUC) to a low-copy one 
(pSC101), so it stands to reason that a different promoter/RBS may be needed in 
this context.  To see if such a promoter/RBS already exists, 

I searched AddGene for a pSC101+AmpR vector, thinking that this problem may 
have laready been solved.  I found pTHSSe_1 (Addgene #109233) from the Voigt 
lab, which does have a different promoter than mine.  This plasmid is also from 
a paper on designing promoters that give constant expression regardless of copy 
number/growth rate/media conditions/etc [SegallShapiro2018]_.  I might want to 
use this idea, if not for every gene, at least for the Zif268 fusion.  I 
imagine that reducing the noise it the expression level of Zif268 could reduce 
the noise of the overall assay.

That said, my real goal is to keep the ratio of Zif268 to target sites 
constant.  That argues for using an unbuffered promoter, as variations in copy 
number would be met by ≈proportional variations in protein expression.  
(According to [SegallShapiro2018]_, expression is proportional to copy number.)  
I could also use an internal control to account for changes in copy number in 
general.

Predicted promoter strengths
============================
I used the Salis Lab "Promoter Calculator" to compare the strengths of the 
promoters in p186 (terminator before AmpR), p189 (no terminator before AmpR), 
and pTHSSe_1 (pSC101 backbone from Addgene with AmpR).  Note that there are 
currently two version of the Salis Lab Promoter Calculator, one that I can't 
easily find a reference for and one that uses a free energy/ML model 
[LaFleur2021]_.  For my plasmids, both versions seem to pretty much return the 
same results, although the values have different scales.  I'll use the free 
energy/ML model, because it seems newer and I think the approach makes sense.

The predictions for the free energy/ML model have the mathematical form of 
:math:`log(\frac{T_X}{T_{"x,ref"}}`, where :math:`T_X` is a transcription rate 
and :math:`T_{"x,ref"}` is the transcription rate for an reference promoter.  
The reference rate is unknown, but I can still sum the rates of promoters 
facing the same direction.  The result will still be relative to the unknown 
reference rate.

- The AmpR promoter (in p186 and p189) has a predicted strength of 312 au, 
  which I think is reasonably strong.

- p186 has a terminator immediately upstream of the AmpR promoter, so 
  presumably that promoter is the only one responsible for expression.

- p189 has 3 additional promoters facing towards AmpR before the nearest 
  terminator, with strengths of 326 (lacUVm), 382, and 324 au.  All of these 
  are at least slightly stronger than the AmpR promoter, so I think it's 
  reasonable to think that AmpR expression would be 4-5x higher with this 
  plasmid.

- pTHSSe_1 has 2 promoters facing towards AmpR before the nearest terminator:

  - P876, 301 au: This is a fragment of the AmpR promoter/RBS from pUC vectors.  
    It's weaker than the full promoter (312 au), but not by much.  The 
    predicted TSS is just after the start codon, though, so I'm not sure if 
    these transcripts will actually be expressed.  It's possible that this 
    fragment just serves as the RBS for the P998 transcript.

  - P998, 419 au: This promoter falls in the tnpA fragment identified by 
    PLannotate.  I really do wonder where that fragment came from, because I 
    expect it's important.

- If I sum the rates for all the relevant promoters for each construct 
  (assuming that the promoters don't interfere with each other), I see:

  ========  =========  =========
  Promoter  Rate (au)  Rel. Rate
  ========  =========  =========
  p186            312       1.00
  p189           1344       4.31
  pTHSSe_1        720       2.31
  ========  =========  =========

  So the pTHSSe_1 promoters by themselves may not fully replicate p189.

Predicted RBS strengths
=======================
Another way to increase expression is to improve the RBS.  I used the Salis Lab 
RBS calculator tool to predict the strength of the existing RBS:

- Sequence: AAAUGCUUCAAUAAUAUUGAAAAAGGAAGAGUAUGAGUAUUCAACAUUUCCGUGUCGCCCUUAUUCC

  - Note that this includes a good amount of the CDS.

- Strength: 1079.30 au
  
For comparison, that's the strongest RBS on the transcript, although there's an 
internal RBS (which appears about 65 aa into AmpR) with a strength of 438.71 
au.

Promoter/RBS designs
====================
Considering the copy numbers and the predicted transcription rates, I probably 
want to increase AmpR expression between 4-100x.  10x might be a reasonable 
compromise, if I need to pick just one number.

I also want to include some sequences straight from pTHSSe_1.

- p186 with tpnA and AmpR promoters from pTHSSe_1

  - Through the predicted beginning of P1002.

- p186 with tpnA promoter from pTHSSe_1 and full AmpR promoter

  - Don't delete anything relative to the above design.  Just insert the 
    missing AmpR promoter sequence.

- p186 with ≈1344 au promoter

  - De Novo DNA was able to design several ≈1344 au promoters, but all were 
    predicted to be ≈100x stronger in the reverse direction.  I don't trust 
    these results.

  - I browsed the parts toolbox for promoters with strengths of about ≈1000 au 
    (although I'm not totally sure that this is the same scale).  pSH003 (A 
    derivative of J23115) has a strength of 1346.10 au; almost exactly what I'm 
    looking for.  Here's the sequence::

      TTTATACGGTTCTTACGAAATAATACAATGGCTTTA

    Once I decide how to attach this to an RBS, I'll plug it back into the 
    promoter calculator to see how strong it's expected to be.

  - Use entire RBS as identified by RBS calculator.  Insert the promoter 
    immediately afterwards.

- p186 with 10-100x RBS.

  - I used the Salis Lab RBS Calculator to design RBSs with a target expression 
    rates of 10x (the single number mentioned above) and 100x (to match pUC) 
    the existing RBS:

    .. datatable:: rbs_designs.xlsx

    It's possible that using these RBSs will affect how well the promoter 
    works, since the initial transcribed region (ITR) is part of the promoter 
    model.

  - The existing RBS (predicted by the Salis Lab RBS calculator) begins exactly 
    at the TSS of the existing promoter (predicted by the Salis Lab promoter 
    calculator).  That's convenient, because it means I can just replace the 
    existing RBS with the designed ones.

Cloning
=======
AmpR is flanked by SR086 (which is the reverse complement of SR071!) in all of 
my constructs.  So it will be easy to use that for Gibson assemblies.
