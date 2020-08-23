************
Optimize RBS
************

I recently realized that I had designed my constructs using a eukaryotic RBS, 
which of course is not compatible with the bacterial IVTT systems (e.g.  
PURExpress) that I'm using.  However, I had a hard time choosing an appropriate 
RBS to use, for the following reasons:

- I encountered conflicting advice on how to choose an RBS.

- RBS strength is known to be context dependent.
  
- I'd prefer to use standard parts, to make my cloning easier to understand to 
  an outsider.
  
- I'd prefer to use the shortest possible sequence, to limit the amount of DNA 
  in my assays.

Ultimately, I'll have to test a number of RBSs if I want to find one that best 
satisfies all these criteria.
  
Background
==========
The most important part of a bacterial RBS is the Shine-Dalgarno (SD) sequence.  
This sequence binds to part of the 16S rRNA with the sequence ``ACCUCCUUA`` 
[Shine1974]_.  A perfectly complementary SD sequence would therefore be 
``TAAGGAGGT``, but most SD sequences only match some of these positions 
(usually matching better in the middle than the edges).  Ideally, the SD 
sequence would end 5 bp before the start codon, e.g.  ``TAAGGAGGTnnnnnATG``, 
but other separations are tolerated.

Considerations
==============
The PURExpress manual makes the following recommendations regarding the RBS:

- The DHFR control plasmid uses the g10-L RBS [Olins1989]_, i.e. the "leader" 
  (L) sequence of gene 10 (g10) of the T7 bacteriophage.  This is a natural 
  sequence.  The g10-L SD sequence is ``tAAGGAGaT`` (with mismatches in lower 
  case).  [Olins1989]_ show that g10-L includes an additional 9 bp upstream of 
  the SD sequence --- ``TTAACTTTA`` --- that are complementary with a 
  completely different region of the 16S rRNA.  This 9 bp sequence is given the 
  name "epsilon", and contributes to strong translation initiation.

  [Olins1989]_ show that a hairpin which forms upstream of the epsilon and SD 
  sequences is not important for translation initiation.  Accordingly, this 
  stem is not included in the DHFR control plasmid.  However, 17 bp from the 
  wildtype g10-L (including an XbaI site) upstream of the epsilon sequence are 
  retained, despite having no clear function.

- NEB recommends the following 5'-UTR for plasmid-based genes: ``[>5 N] AAGGAG 
  [5-8 N] ATG``.  There are a few things that make me skeptical of this 
  recommendation:
  
   - ``AAGGAG`` is (part of) the SD sequence from the DHFR control plasmid, but 
     is not a canonical SD sequence.  This indicates that this recommendation 
     wasn't really optimized, because it's unlikely that the imperfect, 
     incomplete, natural sequence from the control plasmid is the optimal RBS.

   - Similarly, it's unlikely that any base pair before or after the ``AAGGAG`` 
     SD sequence would work equally well, as NEB implies.  Sequences with 
     additional complementarity to the 16s rRNA should work better.

   - The displacement between the RBS and the start codon is vague.  The ideal 
     separation is 5 bp.  That corresponds to 7 bp here, though, because the SD 
     sequence is 2 bp shorter than it could be.  Still, there's no reason to be 
     vague about this.

  All that said, it may be that none of this matters very much.  But it still 
  bothers me that NEB is making recommendations that don't seem to be very well 
  thought-out.

- NEB recommends the following 5'-UTR genes being amplified by PCR: 
  ``gcgaatTAATACGACTCACTATAgggcttaagtatAAGGAGGaaaaaat`` (the T7 promoter and SD 
  sequence are upper case).  Note that this SD has 7 consecutive bp of 
  complementarity with the 16S rRNA, compared with just 6 for g10-L (although 
  NEB doesn't consider the 3' G to be part of the SD sequence---this seems to 
  be a mistake).
  
  I don't know where the ``cttaagtat`` 5' of the SD sequence came from.  It's 
  different than the sequence in the DHFR control plasmid, and an internet 
  search doesn't produce any obvious origins.  It's also notably longer than 
  the 5 bp minimum suggested by NEB above, which is interesting because this 
  sequence is presumably optimized to be short (as it's already on the 
  long-side for a PCR primer).

  It's also worth noting that the `aaaaaat` 3' of the SD sequence creates a 6 
  bp spacing, 1 bp longer than ideal.  Again, it's not clear what the rationale 
  for this is.

All told, I get the impression that NEB did not put much care into their 
recommendations.  I'm sure the recommended RBSs would work fine, but PURExpress 
is expensive, so it seems worth using a good RBS to get the most possible 
expression out of each reaction.

For this reason, I decided to try using the Salis lab RBS calculator to design 
RBSs that satisfy a thermodynamic model.

  
   
In addition to these recommendations, 

g10L
----
T


g10L

