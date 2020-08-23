*************************************
Purify Zif268-repA via rRNA digestion
*************************************

Given that Zif268-repA bound to ORI should be retained by the 100K spin 
filters, I could use those filters to purify and concentrate that complex if I 
can remove everything else bigger than 100 kDa from the reaction.  This 
primarily means removing the ribosomes.  There are several other components of 
the PURE reaction that are nearly 100 kDa [Shimizu2001]_, including notably T7 
RNA polymerase (99 kDa).  However, these components are His-tagged (and were 
purified in the first place using Ni-NTA), so they can be removed using Ni-NTA 
resin.

.. datatable:: t7_rnap.xlsx

One way to get rid of the ribosomes is to digest the rRNA so that the 
r-proteins dissociate and fit through the 100K filter.  All of the rproteins 
individually are small enough to fit through the filter:

.. note::
   
   Another way to accomplish a similar goal is to just dissociate the ribosome.  
   See [Arai2020]_, who (I think) do this in the context of cDNA display by 
   first incubating with 900 mM KCl and 75 mM MgCl₂, followed by incubating 
   with 45 mM EDTA.  Writing this out now, I actually don't see the point of 
   adding EDTA, if the amount of edta isn't even going to exceed the amount of 
   salt added in the previous step.

   On a related note, [Miall1969]_ has evidence that the ribosome begins to 
   dissociates with EDTA concentrations over 2 mM.  They also mention that 
   dialysis into buffer lacking magnesium can dissociate the ribosome.

.. datatable:: rproteins.xlsx

   Sizes of all the E. coli ribosomal proteins listed on `Wikipedia 
   <https://en.wikipedia.org/wiki/Ribosomal_protein#Table_of_ribosomal_proteins>`.

It's not clear which RNases would be most successful at completely digesting 
the rRNA, much of which is structured and buried.  I'm going to start with a 
cocktail of RNase A and RNase T1 (Invitrogen AM2286).  Both are common enzymes, 
and should be capable of digesting RNA nearly to nucleotides.  They are also 
both small enough to fit through the 100K filter:

.. datatable:: rnases.xlsx

Fitzy mentioned that RNase III, which cuts dsRNA and is involved in rRNA 
maturation, might also be useful.  RNase V1 also digests dsRNA, but 
unfortunately seems to have been discontinued by all the manufacturers.

This protocol would also digest the mRNA present in the sample, which is a 
bonus.

Methods
=======

2019/06/28
----------
.. protocol:: 20190628_purexpress.txt

   - DHFR wouldn't be a great control, because the protocol depends on 
     capturing large proteins in the spin filter.  But it would at least tell 
     me that the reaction is basically functioning.

   - The lanes I'd run:

      - crude  (forgot to save an aliquot for this one...ugh)
      - digest
      - 100K flow-through
      - 100K retentate

   - Constructs I'd run:

      - DHFR
      - 11

.. protocol::

   1. Decide how to use RNase cocktail:

      - Invitrogen recommends 2.5 μL for a 50 μL miniprep.
         
         - I can probably use less, becuase I don't have as much RNA.

         - I'm going to use 0.5 μL/rxn, which for my 10 μL reactions is the 
           same proportion as above.

      - No specific buffer recommended or provided.

         - Storage buffer is:

            - 10mM HEPES, pH 7.2
            - 20 mM NaCl
            - 0.1% Triton
            - 1 mM EDTA
         
         - For comparison, 1x PBS is:

            - 10 mM NaH₂PO₄, ph 7.4
            - 1.8 mM KH₂PO₄, ph 7.4
            - 137 mM NaCl
            - 2.7 mM KCl

         - I couldn't find the composition of NEB's PURExpress buffer, but the 
           buffer used in the original PURE system [Shimizu2001]_ is (50 μL for 
           mass/quantity units):

            - 9 mM magnesium acetate
            - 5 mM potassium phosphate, pH 7.3
            - 95 mM potassium glutamate
            - 5 mM ammonium chloride
            - 0.5 mM calcium chloride
            - 1 mM spermidine
            - 8 mM putrescine
            - 1 mM dithiothreitol (DTT)
            - 2 mM each ATP and GTP
            - 1 mM each of CTP and UTP
            - 10 mM creatine phosphate
            - 2.8 A260 units tRNA mix (Roche, Mannheim, Germany)
            - 0.5 μg 10-formyl-5,6,7,8-tetrahydrofolic acid
            - 0.1 mM each of amino acids
            - 12 pmol (32.4 μg) ribosome
            - 1 μg IF1
            - 2 μg IF2
            - 0.75 μg IF3
            - 1 μg EF-G
            - 2 μg EF-Tu
            - 1 μg EF-Ts
            - 0.5 μg RF1
            - 0.5 μg RF3
            - 0.5 μg RRF
            - 30–300 units of each ARS and MTF
            - 0.2 μg creatine kinase (Roche)
            - 0.15 μg myokinase (Sigma, St. Louis, MO)
            - 0.054 μg nucleoside-diphosphate kinase
            - 0.1 units pyrophosphatase (Sigma)
            - 0.5 μg T7 RNA polymerase

         - I'm just going to add RNase directly to reaction, that's basically 
           what's recommended anyways.

         - Theres about 10 μg of non-ribosome protein components, not counting 
           the ARSs and MTFs, because I'm not sure what a unit is.

      - No time or temperature recommendations made.

         - Invitrogen says the enzymes can be added to a "normal" restriction 
           digest.
           
         - I took that to mean that 37°C for 15 min would be reasonable.

         - If it doesn't work well, I might try longer times (or more enzyme).

   2. Add 0.5 μL RNase cocktail to each reaction.

   3. Dilute reactions to 100 μL with 1x PBS.

   4. Take 10 μL aliquot.

   5. Load 100K spin filters.

   6. Spin 30 min, 1500g, 4°C

   7. Save 10 μL flow-through (ft) and retentate (ret) aliquots.

   8. Repeat steps 3-7.

   9. Repeat steps 3-7 again.

   10. Run SDS-PAGE.

.. figure:: 20190701_rnase_digestion.svg

- I forgot to save an aliquot of the crude reaction, but comparing to the crude 
  reaction from another gel (`20190626_purify_controls.svg`), I don't see any 
  differences.  I don't even see the RNases, although at ~10 and 13 kDa, they 
  might have run off the bottom.

- The DHFR control is not retained by the filter.  Some more comes off on the 
  second wash, but this is probably the protein that was in the dead volume of 
  the spin filter the first time.

- Zif268-repA (11) is retained by the filter for all 3 spins, with no 
  detectable protein passing through the filter.  This suggests that repA is 
  bound to DNA as it should be, otherwise the Zif268-repA protein should fit 
  through the filter.  It also means that I can concentrate the complex, which 
  will make things easier.

- The ribosome seems to be mostly retained by the filter, suggesting that the 
  RNase treatment was not sufficient to destroy the ribosomes.  Extra spins do 
  not seem to get rid of any more of the ribosomes.

- I added RNase inhibitor to the IVTT reaction, which could obviously interfere 
  with the RNase digestion.  RNase inhibitor is necessary when using 
  miniprepped template DNA (since RNase is added during a miniprep), but 
  shouldn't be as necessary with PCR-amplified template.

Results
=======
I do not think RNase digestion is a promising way to get rid of the ribosomes, 
but if I were to look into this again:
  
- Leave out the RNase inhibitor.
- Use more RNase and for a longer time.  At least 1h, maybe overnight.
- I don't think I need any more than 2 spins.

