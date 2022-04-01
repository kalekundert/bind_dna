****************************************
Remove unannotated/unnecessary sequences
****************************************

2022/03/09:

I have been able to successfully assemble the combined B1H plasmid, but the 
approach that ended up being successful left a lot of sequences that are either 
unannotated or annotated and expected to be unimportant.  For the sake of 
simplifying things before I do a lot of work making libraries, I want to see if 
I can remove any of these sequences.

Here are the sequences I'm thinking of removing:

- f1 ori

  - I believe this is only used to make phage.

- unannotated sequence between Rep101 and f1 ori

- unannotated sequence between AmpR and pSC101 ori

- GAL1 sequence 

  - This is a (fragment of a) yeast coding sequence, presumably related to the 
    (in yeast) inducible GAL promoter.

Notes:

- I'm going to use p194 as a template.  p193 is only different in that the f1 
  ori is split up.  This difference won't matter if I can delete the f1 ori, 
  and having the f1 ori all in one place might make the results more clear.

- I'm going to ignore p186 for now.

- When designing primers to delete the GAL11 fragment, I ran into problems with 
  some problems with the 5' UTR of the HIS3/URA3 transcript having the same 
  sequence.  If I decide to use a [Mutalik2013]_ promoter/RBS for the rpoZ 
  fusion, I can probably eliminate this duplication.

2022/03/16:

I was able to get colonies for all of the deletions except p203, which deletes 
the unannotated sequence between AmpR and the pSC101 ORI.  I repeated the PCR 
with a temperature gradient, and found very good amplification at a range of 
temperatures.  

.. protocol:: 20220315_pcr.txt

.. figure:: 20220315_optimize_ta_p203.svg

I didn't have space on the gel to include the original PCR I ran, but it was 
done with :math:`T_A = \pu{64°C}`, so there should've been some product.  I 
repeated the cloning with the 58°C reaction, and once again got no colonies.  
Given that the PCR product looks perfect, and that I've successfully used my 
KLD mix a lot recently, I think there may be some essential sequence in this 
region.

That said, according to [Sugiura1993]_, the annotated ORI includes all 
essential sequence.  In fact, they deleted the ORI at the same exact place as I 
did (Δ3) and saw no decrease in fitness.

Here are the options I can think of:

- Nothing's wrong with the plasmid, my cloning just isn't working for 
  mysterious reasons.

- I need more sequence upstream of the ORI.  This is contrary to what 
  [Sugiura1993]_ reported, but maybe there's something pathological about the 
  upstream AmpR sequence that requires more buffer.

- I need more sequence downstream of AmpR.  This is a little hard to believe, 
  because I have a stop codon and a terminator already, plus 20 bp of 
  downstream context.  There is a somewhat terminator-looking sequence 195-222 
  bp downstream of the λ t0 terminator, but it seems unlikely that it's 
  essential.


