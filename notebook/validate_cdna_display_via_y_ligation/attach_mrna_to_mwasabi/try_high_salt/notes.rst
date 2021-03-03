*************
Try high salt
*************
Most mRNA protocols (with the strange exception of [Barendt2013]_) call for a 
high-salt incubation after protein expression, and claim that this step is 
essential for protein/mRNA coupling.  For example, [Liu2000]_ states:

  By the end of a normal translation reaction, only a small amount of the in 
  vitro-synthesized protein has been converted to its mRNA-protein fusion [12].  
  [...]  Subsequent experiments clearly demonstrate that three factors can 
  routinely provide for more than 40% of the input RNA and 50-80% of in vitro- 
  translated protein to be synthesized as fusions: (1) posttranslational 
  addition of Mg2+ and K+, (2) use of a flexible linker of the correct length, 
  and (3) incubation after completion of translation at room temperature if K+ 
  and Mg2+ are added or long incubations (12-48 hr) at low temperature (-20°C) 
  if they are not.
      
The precise parameters of the incubation (salts, concentrations, times) need to 
be optimized for each protein, and therefore vary somewhat between all the 
protocols.  Here I will test if such an incubation will allow me to observe 
coupling.

Considerations
==============
Reagents:

- 3M KCl
- 1M MgOAc

I'll try using the salt concentrations from [Naimudden2016]_ first.  All of the 
references are pretty similar, but [Naimudden2016]_ has the highest salt 
concentrations.  Based on the plots in [Liu2000]_ (pg 282), it's better to err 
on the side of having too much salt.  Plus, [Naimudden2016]_ is the method I'm 
most trying the replicate.

.. note::

  I would've followed the protocol from [Barendt2013]_, since that's the only 
  method with the same ribosome as I'm using, but it doesn't mention a 
  high-salt incubation.  I think this might be an omission.

[Reyes2021]_
------------
.. note::

  I didn't find this paper until 2021/02/26.

- Expression reaction:

  - PUREfrex 1.0

    - "According to GeneFrontier, PUREfrex 2.0 is optimized for maximizing the 
      protein translation efficiency and therefore it does not match the 
      purpose of stabilizing the ribosome at the mRNA/DNA‐tag junction to form 
      an mRNA‐peptide conjugate."

    - I have PUREfrex 2.0, but maybe it's worth trying PUREfrex 1.0 as well.

  - 37°C for 30 min

- Coupling reaction:

  - 32.5 mM MgCl₂
  - 375 mM KCl

- Electrophoresis:

  - 2x Laemmli buffer
  - Spin 1m, 10000g to remove salt precipitate

[Liu2000]_
----------
- Expression reaction:

  - Incubate at 30°C for 30-90 min.
  - Incubate at 4°C for 40 min.

- Coupling reaction:

  - 50-100 mM Mg²⁺
  - 300-600 mM K⁺
  - Incubate at 25°C for 1h.

[Keefe2001]_
------------
- Expression reaction:

  - 200-800 nM mRNA/linker
  - 100 mM KCl
  - 500 µM MgOAc
  - Incubate 1h at 30°C

- Coupling reaction:

  - 50 mM MgCl₂
  - 565 mM KCl
  - Incubate 5 min at 25°C, or up to 1 week at -20°C

[Cotten2011]_
-------------
- Expression reaction:

  - 50-120 mM Kcl
  - 0.3-1.0 mM MgOAc
  - 200 nM mRNA/linker
  - Retic lysate
  - Incubate at 30°C for 90 min

- Coupling reaction:

  - 50 mM MgCl₂
  - 580 mM KCl
  - Incubate at 25°C for 30 min.
  - Incubate at -20°C overnight.

[Barendt2013]_
--------------
- Expression reaction:

  - PURExpress
  - 60 pmol mRNA
  - Incubate at 37°C for 30 min
  - Incubate at 25°C for 10 min

[Naimudden2016]_ --- Wheat germ extract
---------------------------------------
- Incubate mRNA/linker at 65°C for 10 min, then chill on ice

- Expression reaction:

  - Wheat germ extract
  - 50 mM KOAc
  - Incubate at 25°C for 10 min

- Coupling reaction:

  - 65 mM MgCl₂
  - 750 mM KCl
  - Incubate at 25°C for 1h.

[Naimudden2016]_ --- Rabbit reticulocyte lysate
-----------------------------------------------
- Incubate mRNA/linker at 65°C for 10 min, then chill on ice

- Expression reaction:

  - rabbit reticulocyte lysate
  - salt unspecified
  - Incubate at 30°C for 10 min.

- Coupling reaction:

  - 65 mM MgCl₂
  - 750 mM KCl
  - Incubate at 37°C for 2h.



Results
=======

.. protocol:: 20200908_rnasezap_purexpress_couple_puromycin_gel_laser_scanner.txt

.. figure:: 20200908_compare_salt_linker.svg

- As in :expt:`65`, the mRNA seems to have been fully degraded.  
  
  This was a shorter expression reaction (30 min at 37°C), although it was 
  followed by a +/- salt incubation (1h at 25°C with high salt) that I didn't 
  do previously.  Overall, though, this was a shorter incubation at a lower 
  temperature.  But the apparent degradation is unchanged.  I'm interested to 
  do a timecourse, and control +/- each component of PURExpress.

- The GFP signal in this gel was much fainter than usual.  Two possible 
  explanations:

  - Less protein was expressed because the reaction was shorter (30 min instead 
    of 2h).

  - The mWasabi didn't have enough time to mature.  I checked FPbase for the 
    maturation time of mWasabi, but it wasn't listed.

- The KCl from the high-salt incubation precipitated the SDS in the loading 
  buffer.  I ran the gel anyways and it seemed to turn out fine, but it's 
  something to be aware of.  I might try running a urea gel in the future.  

- The combination of only using 0.6x linker and doing the urea washes was very 
  effective at removing unligated linker.
