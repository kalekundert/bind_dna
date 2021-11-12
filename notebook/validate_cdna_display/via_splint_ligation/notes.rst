*******************
Via splint-ligation
*******************

I'm suspicious that RNase H contamination in PURExpress is preventing 
:expt:`53` from working, due to the RNA-DNA duplex inherent to Y-ligation.  One 
possible way to get around this is to replace Y-ligation with splint ligation.

Splint ligation is the technique normally used in mRNA display to attach 
puromycin to the mRNA.  [Naimudden2016]_ used Y-ligation instead for cDNA 
display, claiming that the latter is more efficient.  In order to use splint 
ligation for cDNA display, I'd need to modify the oligos to add an RT primer, 
but this should be doable.

My thought is to begin by trying to replicate [Barendt2013]_, and then to worry 
about adding the RT primer.

Considerations
==============

RT primer
---------
- Sequester RT primer with dsDNA, use exonuclease to free:

  - Modify desired strand with 5' biotin
  - Treat dsDNA with 5'→3' exonuclease that would be blocked by biotin (e.g. λ 
    exonuclease, T7 exonuclease).
  - Can also use biotin in affinity selection steps.

  - Biotin would also block RT, though, right?

DNA vs RNA linker
-----------------
- I'd prefer to have a DNA linker:

  - Cheaper
  - No RNA in final construct; should be more stable.

- DNA *can* be used as the template for MMLV RT:

  https://www.neb.com/faqs/0001/01/01/can-dna-be-used-as-a-template-for-m-mulv-reverse-transcriptase

  - This suggests that the following architecture might work:

    .. figure:: architecture.svg

    The reverse transcriptase would have to transcribe a few bases of DNA 
    before getting the the RNA, but it might work.
