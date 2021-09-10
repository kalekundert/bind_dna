******************
Optimize PCR yield
******************

:expt:`127` will require extremely large amounts of PCR product.  In order to 
make this more achievable, I want to optimize the amount of product I can get 
from a single PCR reaction.

.. protocol:: 20210827_optimize_pcr_yield.pdf 20210827_optimize_pcr_yield.txt

.. datatable:: 20210827_optimize_pcr_yield.xlsx

Observations:

- 35 cycles is not enough to complete the reaction, although additional cycles 
  only modestly improve yield.

- I'd like to see if I can saturate the reaction (e.g. use up all the free 
  nucleotides/ATP/whatever component of the master mix is limiting).  The next 
  step in this direction would probably be to try a reaction with 55 cycles.

  According the `Thermo 
  <https://www.thermofisher.com/us/en/home/life-science/cloning/cloning-learning-center/invitrogen-school-of-molecular-biology/pcr-education/pcr-reagents-enzymes/pcr-cycling-considerations.html>`_:

    More than 45 cycles is not recommended as nonspecific bands start to appear 
    with higher numbers of cycles.

  If I try 55 cycles, it'd be important to run a gel to make sure the product 
  is reasonably high quality.

- None of the reactions even get close to using up all of the available primer.  
  Despite that, greater primer concentrations do give more product for the 45 
  cycle condition.
