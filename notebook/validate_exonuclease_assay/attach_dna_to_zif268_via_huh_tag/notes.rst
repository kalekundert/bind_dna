********************************
Attach DNA to Zif268 via HUH-tag
********************************

See :expt:`47` for the motivation behind using HUH-tags.  Here I want to show 
that I can attach DNA to Zif268.  I don't expect the Zif268-PCV2 fusion to 
behave differently than dCas9-PCV2.

Results
=======

Cloning --- 2020/08/05
----------------------
2020/08/05:

.. protocol:: 20200805_make.txt

2020/08/11:

I'm having trouble amplifying f45.  I'm going to check my primers and try a 
temperature gradient.  I should also sequence p106, since I got it from Jorge 
and I haven't confirmed its sequence myself.

.. update:: 2020/08/17
   
   There was an E→G mutation in the linker peptide that was causing my primers 
   to not bind.  I don't think the mutation is a problem, so I  just ordered 
   new primers to account for the actual sequence.

2020/08/19:

.. protocol:: 20200819_pcr.txt

.. figure:: 20200819_optimize_ta_f45.svg

I still couldn't amplify f45 even with the new primers, so I did a gradient 
PCR.  Interestingly, it seems as if annealing temperatures both higher and 
lower than 66°C (the NEB recommendation, and what I had been using previously) 
would work.  I don't think I'll have to do another PCR, though: I should have 
enough product left over from these ones.

