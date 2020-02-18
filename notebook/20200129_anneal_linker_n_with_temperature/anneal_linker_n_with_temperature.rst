********************************
Anneal linker-N with temperature
********************************

I was able to find a number of different temperature protocols for annealing 
oligos, see :expt:`20191209_ligate_linker_n_via_naimudden2016`.  They're 
similar for the most part---start hot, gradually cool---and I generally get the 
impression that any would be fine, but I want to test to see if there's any 
reason to prefer the longer, slower protocols.

Results
=======
Below are the controls I want to run:

- mRNA (f11) only, hold at 4°C
- linker-N only, hold at 4°C
- no PBS, hold at 4°C

Below are the temperature protocols I want to test (from
:expt:`20191209_ligate_linker_n_via_naimudden2016`):

- Hold at 4°C.
- 95°C for 2 min, leave on bench until cool.
- 95°C for 2 min, 95°C→25°C over 45 min, hold at 4°C
- 95°C for 5 min, 95°C→25°C over 70 min, hold at 4°C
- 95°C for 5 min, 95°C→Tm at 1°C/min, Tm for 30 min, Tm→25°C at 1°C/min, hold 
  at 4°C

The Tm of the Y-tag sequence ``CAAGGGCGGGGGGCGGCGGGG`` is 88°C according to NEB 
TmCalc and 75–76°C according to SnapGene.  I'm inclined to trust SnapGene more, 
because `they claim <https://www.snapgene.com/support/faq/>`_ that their Tm 
algorithm is more modern/accurate than NEB's, plus I'm not using Q5 here.

Based on :expt:`20191219_anneal_linker_n_with_salt`, I will use PBS for all 
reactions.

.. protocol::

   See binder, 2020/01/30.  By mistake I added 0.4 µL too much water to each 
   reaction (4.4 µL total), so everything is a bit too dilute.

