**********
Clone p243
**********

2022/10/05
==========

See :expt:`202`.

2022/10/17
==========
Interpreting the sequencing results for the "f203" fragment I purified above:

- f203 has the exact same sequence as p242.

  - Either the BsmBI digestion of p242 failed, or the two fragments from that 
    digestion religated.  Neither is consistent with the fact that I did a gel 
    purification, and had a ≈900bp band to purify (although the gel was 
    smeary).

    - Would Sanger sequencing of f201 be informative?  It should tell me if the 
      digestion worked, which would distinguish between the two possibilities 
      presented above.

    - Maybe this was another partial digestion.  That could explain why 
      re-ligation happened.

- Next steps:

  - Run PAGE gel of f201 (and p242 digest?)
  - PAGE-purify f201?

  - Sequence f201 and f200?

2022/11/14
==========

.. protocol:: 20221114_make_f200.pdf 20221114_make_f201.pdf 20221114_make_f200.txt 20221114_make_f201.txt

.. figure:: 20221114_digest_p239.svg

- I don't think that p239 has the right sequence.

  - Sanger sequencing on p239 failed, see :expt:`200`.

  - Despite this, I thought that p239 probably had the right sequence, because 
    (i) I purified the two fragments so carefully and (ii) I got >10x more 
    colonies with insert than without. 

  - However, Esp3I seems to have no effect on the plasmid, which suggests that 
    the cloning did not in fact work.

  - f200 should be 3.4 kb, but there is no band of that size.  Some 
    explanations:

    - The plasmid contains no Esp3I sites

      - In this case, the band would either be supercoiled plasmid or (somehow) 
        linear DNA.

    - My Esp3I is completely inactive.

      - This would have the same symptoms as above, but is unlikely because I 
        successfully used this enzyme just last week.

    - There's a lot of smear.  That could be RNA and/or genomic DNA.

.. figure:: 20221114_make_f201.svg

- No yield from f210 PCR reaction.

2022/11/15
==========
The f210 PCR reaction didn't produce any product yesterday, so today I did a 
:math:`T_A` optimization:

.. protocol:: 20221115_optimize_tm_f210.txt

.. figure:: 20221115_optimize_ta_f210.svg

.. datatable:: 20221115_optimize_ta_f210.xlsx

2022/11/16
==========
.. protocol:: 20221116_make_f200.pdf 20221116_make_f200.txt 20221116_make_f201.pdf 20221116_make_f201.txt

- f201 digestion did not go to completion.

  - I only used a 3.33x excess of enzyme, instead of the usual 10x excess.  I 
    might need to do 10x.

- f200 digestion seemed to have two bands, very closely spaced.

  - Not sure what this could be.
    
  - Tried to mostly cut out the upper band, since it was closer to 3.4 kb.  Not 
    sure this was the right decision, though.  Maybe I should've just gotten 
    both.

  - Maybe next time I should just do an overnight digestion, and no gel 
    purification.  I didn't plan to do a gel purification, but I my E-Gel made 
    it look like there was uncut plasmid I had to remove.  The purification gel 
    didn't look so much like that, though.

2022/11/18
==========
.. datatable:: 20221118_electrotransform_p243.csv

- Got 25x more transformants with insert than without, implying that 96% of the 
  transformants have the insert.
  
- Sanger sequencing failed, despite the fact that I measured the concentration 
  of the plasmid by Qubit.




