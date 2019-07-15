************************
Optimize qPCR conditions
************************

For my initial control experiments, I'm using qPCR to quantify relative DNA 
concentrations.  Before I can do this, I need to know conditions where the DNA 
in question amplifies well.  The conditions I want to optimize are:

- Annealing temperature (Ta)
- Primer concentration

Annealing temperature (Ta)
==========================
.. protocol:: 20190610_pcr.txt

   - I used SsoAdvanced™ Universal SYBR® Green Supermix (Biorad 1725270) rather 
     than Q5 master mix.

   - I used the two-step thermocycler protocol recommended by BioRad for use 
     with the above polymerase:

      - 95°C for 30s

      - Repeat 40x:

         - 95°C for 10s
         - 55-63°C for 15s
      
      - Melt curve

         - 65-95°C in 5s steps 0.5°C

   - I used a BioRad CFX96 thermocycler with a white plate and a 
     pressure-activated optical seal.

.. figure:: 20190610_optimize_ta.svg

   Cq values for different annealing temperatures.

- Ta=55.5°C gave the lowest Cq values in this experiment.  However, all Ta<59°C 
  gave pretty similar values, and there wasn't a clear minimum.  I might be 
  inclined to use 58°C in the future, since it's closer to the 60°C "optimal" 
  temperature for SsoAdvanced.

- The melt curves are not shown here, but they were all unimodal and gave no 
  sign of primer dimers or nonspecific amplification.  This isn't too 
  surprising, because the template itself is pretty clean, and I've already 
  seen pretty clean amplification with these primers in regular PCR.

Primer concentration
====================
.. protocol::

   See binder.

.. figure:: 20190610_optimize_primer_conc.svg

   Heatmap showing the Cq values for the different primer concentrations.  The 
   standard deviations are not show, and are large in some cases (up to 1.62).  
   But the standard deviation for the 500/500 condition is low (0.16).
   
- There isn't much difference between the various primer concentrations, but 
  500 nM of both primers gives the lowest Cq.  It makes some sense that the 
  highest primer concentration would work the best with very pure template.  
  I'll probably continue using 500 nM for both primers moving forward.

- It's interesting that the Cqs in this experiment are so much higher (2 Cq 
  units) than the previous experiment.  These were done on the same day, with 
  the same dilutions of the same reagents.  The 500 nM, 500 nM condition is an 
  exact replicate from the first experiment.  It makes me wonder if I made some 
  sort of mistake...
   
Efficiency
==========
.. protocol:: 20190611_pcr.txt 

   - Serial dilution of template DNA:

      - Initial: 100 pg/µL
      - 7 steps:
         
         - 45 µL water
         - 5 µL previous


