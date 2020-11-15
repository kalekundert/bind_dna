**************
Try longer ori
**************

The primer I've been using to install the ori sequence (o102) has the minimal 
sequence needed to be recognized by PCV2-Rep.  There's plenty of evidence that 
this sequence should work (:expt:`46`), but Jorge also recommended adding 15 
arbitrary nucleotides between the ori sequence and the spacer and I want to 
give that a try.

Considerations
==============

Spacer sequence
---------------
Jorge didn't recommend a specific spacer sequence to use, he just said 
"NNNNNNNNNNNNNNN".  So I do need to pick a sequence.  Probably any choice is 
fine, but here are some options:

- Poly-A
- Random
- Random but controlling for GC content
- Corresponding sequence from PCV2 genome

  I found 4 PCV2 genomes after a very brief search.  All have the same ori 
  sequence and context.  I used the 12-nt ori sequence from [VegaRocha2007]_ 
  and included 6 and 15 nt of upstream context and downstream context, 
  respectively.

  - PCV2 refseq: cgctgtAAGTATTACCAGcgcacttcggcagcg
  - PCV2a:       cgctgtAAGTATTACCAGcgcacttcggcagcg
  - PCV2b:       cgctgtAAGTATTACCAGcgcacttcggcagcg
  - PCV2c:       cgctgtAAGTATTACCAGcgcacttcggcagcg

  After skimming some literature of PCV2 and looking at figure 7 from 
  [VegaRocha2007]_, I realized that it makes sense to include the palindromic 
  stem flanking the ori and the 2 hexamer repeats downstream of that stem.  
  Both features along with the ssDNA hairpin, are implicated in PCV2 binding:

  - gaagtgcgctgtaagtattaccagcgcacttcggcagcggcag

  It's also worth noting that the ori sequence used by Jorge is actually wrong: 
  the real ori begins with "tAAG" not "taAAG".

I'm going to try two sequences: 

- o204: the native genome context
- o205: exactly what Jorge told me (with N replaced by A).


Results
=======

First try --- 2020/11/05
------------------------
.. protocol:: 20201105_compare_ori_linkers.txt

.. figure:: 20201105_compare_ori_linkers.svg

- The highest MW bands (≈10kb) have both protein and DNA, but appear green 
  because the protein signal is very faint.

- The DNA is cleaner in this experiment than it has been previously, e.g.  
  :expt:`46` and :expt:`68`.  This is because I ordered new primers for this 
  experiment, specifically to get cleaner amplification.

- I don't know why the DNA with the A15 ori shifts dramatically in the presence 
  of dCas9-Rep without EDTA.

  - I've seen this same effect in previous experiments, both of which used ≈10x 
    less protein:

    - :expt:`46`: Nearly all DNA in shifted band.
    - :expt:`68`: Only about 20% of the total DNA is in the shifted band.

  The easiest conclusion to make is that the DNA is reacting with truncated 
  protein.  With the EDTA and −protein controls, it's hard to argue that Rep is 
  not involved in the shift.  And there is a *very* faint red band under the 
  shifted green band.

- I don't know what the ≈200 bp DNA species in the A15 lanes is.

  It's slightly bigger than f105, but I wouldn't rule out the possibility of 
  contamination.  This band didn't seem to be present in the E-gel I ran after 
  the PCR (not shown), so it probably isn't an incorrectly-amplified PCR 
  product.

Repeat A15 linker --- 2020/11/12
--------------------------------
.. protocol:: 20201112_repeat_a15_linker.txt

.. figure:: 20201112_repeat_a15_linker.svg

Observations:

- I used a fresh aliquot of protein in this reaction.

- Compared to the previous experiment:

  - The coupled DNA is mostly in the ≈10kb band rather than the ≈350 bp band.

  - There's much less ≈200 bp contaminant (although still some).

Conclusions:

- I think I somehow contaminated the previous (2020/11/05) f107 reactions with 
  f105.  That explains why there's a big ≈200 bp band, and why it would go away 
  when I setup the reaction with only a single template.

  - There is still a very faint band at 200 bp.  That may represent a low level 
    of contamination in my stock.

- I think that all the DNA is reacting with protein, but that (i) the protein 
  is partially degraded by either proteases or free-thaw cycles and (ii) the 
  Rep domain is poorly stained by Coomassie.

  - Clearly the DNA is reacting with something.  Given that EDTA completely 
    eliminates the reaction, it's very likely that that something is Rep.

  - The ≈350 bp band is too small to be the dCas9-Rep fusion, so it must be the 
    free Rep domain (bound to DNA).

  - The only way that free Rep domain could be present is if the dCas9-Rep 
    fusion is being degraded somehow.

    - This is corroborated by the fact that this experiment, which used a fresh 
      aliquot, has much less DNA in the ≈350 bp band.

    - I'm not sure if proteases or freeze-thaw cycles are the more likely 
      culprit.

    - It also follows that some of the unreacted Cas9 I've seen in other 
      experiments had simply lost the Rep domain, and therefore couldn't react.

  - In the previous experiment (see inset), there's a very faint band in the 
    Coomassie channel corresponding to ≈350 bp band in the GelGreen channel.

    - I had originally assumed that this was cross-talk from the GelGreen 
      channel, but that's not the case.  The DNA-only control lanes have no 
      such cross-talk at all.  (Note that there is some cross-talk in the 
      reverse direction, e.g. some signal from the Coomassie channel does bleed 
      into the GelGreen channel.)

    - This band must genuinely represent the Rep domain, then.  The faintness 
      of the band must mean that Coomassie does not bind well to this domain.

      - Coomassie mostly binds to Arg and (to a lesser extent) Tyr and His 
        [Compton1985]_.

      - Proteins typically have only a few strong Coomassie binding sites, and 
        the number of sites is not strongly correlated with the size of the 
        protein [Congdon1993]_.  It may be that a strong binding site requires 
        an arginine and an aromatic residue in close proximity [Congdon1993]_.

      - The PCV2 Rep domain is 248 aa (27.9 kDa) and has 8 Arg, 3 Tyr, and 3 
        His.

.. todo::

  - Try a different staining technique:

    - SYPRO Orange: Binds SDS coat.

    - No-Stain: Covalently reacts with lysine residues.

    - Coomassie R-250: I've been using SimplyBlue SafeStain, which is a G-250 
      (i.e. colloidal) stain.  I've seen some comments suggesting that these 
      two dyes interact differently with proteins and that R-250 is more 
      sensitive, but I can find any good sources to back up these comments.

  - Express protein using PURExpress.

