*************************
Pick initial model system
*************************

I need to pick a DNA binding domain to use for my initial experiments.

Kamesh:

- Zif268
   - a.k.a. Egr1
   - 90 residues
   - 1aay
   - monomeric

Methods
=======

Reverse translate Zif268
------------------------
I found the wildtype Zif268 sequence on pB1H1-Zif268 from Addgene.  However, 
that sequence contains restriction sites that I might want to use, so I reverse 
translated the protein sequence myself using E. coli codons and avoiding those 
restrictions sites::

   $ ./reverse_translate.sh zif268.prot

Control targets
---------------
The cognate sequence for Zif268 is ``GCGTGGGCG``

As a negative control, I decided to use ``GCGAAAGCG``.  This is the negative 
control used by [Bulyk2001]_.

Primer design
-------------
I'm going to clone into a pUC19 vector.  My plan is to do Golden Gate cloning 
with BsmBI and SapI.  This will remove most of the golden gate sites in the 
backbone, which may be useful in the future (one BsaI site will remain).  It 
will also remove 868 bp of nonessential sequence (about 30% of the plasmid), 
which should give me yields that much higher.

For the perspective of my inserts, the sticky ends will be:

- 5': ``5'-CGCG-3'`` (BsmBI)
- 3': ``3'-CGA-5'`` (SapI)

