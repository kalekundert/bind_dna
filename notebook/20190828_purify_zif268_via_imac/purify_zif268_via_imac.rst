**********************
Purify Zif268 via IMAC
**********************

I noticed in :expt:`20190723_confirm_zif268_emsa` that Zif268 might bind 
Ni-NTA.  This observation is supported by literature which suggests that the 
C2H2 motifs in some zinc-fingers can bind Ni-NTA or Zn-NTA [Vorackova2011]_.  
It might be convenient to be able to purify untagged zinc-fingers using Ni-NTA, 
so I want to test more carefully how well this works.

Considerations
==============

[Vorackova2011]_ protocol
-------------------------
.. protocol:: 

   - Resuspend cells:
     
      - 50 mM Tris, pH 8.0
      - 150 mM NaCl
      - 1 mM EDTA
      - 1 mg/mL lysozyme
      - 0.05% β-mercaptoethanol
      - 100 μg/mL PMSF
      - Protease inhibitor

   - Sonicate

   - Add 0.1% deoxycholate

   - Pellet; keep supernatant

   - Wash with:

      - 50 mM Tris, pH 8.0
      - 1000 mM NaCl
      - 1 mM EDTA
      - 0.5% Triton X-100

   - Pellet; keep supernatant

   - Wash with:

      - 50 mM Tris, pH 8.0
      - 1500 mM NaCl
      - 1 mM EDTA
      - 0.5% Triton X-100
      - 0.05% β-mercaptoethanol

   - Pellet; keep supernatant

   - Dialyze supernatants containing protein overnight at 4°C:

      - 50 mM NaPO₄, pH 7.5
      - 500 mM NaCl

   - FPLC: 

      - Charge IMAC columns with either Ni or Zn, following manufacturer's 
        instructions.

      - Load dialyzed protein onto column.

      - Wash with 50 mL dialysis buffer

      - Elute in dialysis buffer + 2M NH₄Cl

   - Dialyze eluted protein:

      - 50 mM NaPO₄, pH 7.5
      - 500 mM NaCl
      - 0.01% β-mercaptoethanol
      - 1 μM ZnCl₂

The elution conditions for this protocol are very harsh, but I don't see why I 
couldn't try eluting with imidazole (or even ZnOAc).

How many beads?
---------------
- PURE reactions (not PURExpress, see 
  :expt:`20190626_purify_zif268_repa_via_rrna_digestion`) have at least 10 μg 
  of His-tagged components in a 50 μL reaction (2 μg in a 10 μL reaction). 
  
   - Note that the amount of some components is given in enzymatic "units", 
     and I don't know how those units correlate to μg.

- PURExpress has 10x more ribosomes than PURE (see 
  :expt:`20190626_purify_zif268_repa_via_rrna_digestion`).  If I assume that 
  PURExpress has 10x more of all the His-tagged components as well, that would 
  be 20 μg of protein in each 10 μL reaction.

- In the PURExpress manual, the example yield calculation gives a value of 
  4.67 μg DHFR for a 25 μL reaction.  If I assume this is representative, I 
  can expect about 2 μg of protein in a 10 μL reaction.

- Putting everything together, I would expect there to be between 4-14 μg of 
  His-tagged protein in a 10 μL reaction.

- That corresponds to 20-50 μL of beads (see: Qiagen magnetic bead manual, 
  page 17).

- Magnetic beads are better than agarose beads for micro-scale 
  purifications.


Methods
=======

Ni-NTA purification --- 2019/09/09
----------------------------------
.. protocol:: 20190830_purexpress.txt

   Setup two reactions: one with no template and the other with 34 as the 
   template.

   ***

   See Qiagen handbook for buffer recipes.  See binder for buffer protocols.

   - Ni-NTA purification

      - Dilute reaction to 250 μL in lysis buffer
      - Save 10 μL aliquots ("crude")
      - All following steps only apply to the 34 reaction (i.e. not the minus 
        template control).
      - Add 50 μL bead solution.
      - Mix continuously for 1h at 4°C
      - Separate beads
      - Save a 10 μL aliquot ("flow-thru")
      - Wash with 500 μL wash buffer
      - Save a 10 μL aliquot ("wash 1")
      - Repeat wash
      - Save a 10 μL aliquot ("wash 2") [KBK: I forgot to do this step]
      - Add 25 μL elution buffer
      - Incubate of ice for 2 min.
      - Remove beads.
      - Save a 10 μL aliquot ("eluate")

   - SDS-PAGE

      - Samples:

         - 10 μL IVTT
         - 3.85 μL 4x buffer
         - 1.54 μL 10 reducing agent
         - 70°C for 10 min

      - 165V, 42 min.


.. figure:: 20190910_purify_zif_wo_his6.svg

- I forgot that the PURExpress reactions contain a bunch of His-tagged 
  components, which are obviously purified by Ni-NTA.  These components 
  actually serve as an internal positive control for this experiment, which is 
  nice.

- I don't think any Zif268 was retained by the Ni-NTA.  A band corresponding to 
  Zif268 is visible in the crude reaction and the flow-thru aliquot.  There 
  appears to be no Zif268 band in the eluate, but it's hard to say for sure 
  because I didn't do the purification on the negative control reaction as a 
  comparison.

  Perhaps this is because I had Zn in the reaction buffer.  In any case, it's 
  probably for the best if my protein doesn't have any weird interactions with 
  Ni-NTA on its own.

