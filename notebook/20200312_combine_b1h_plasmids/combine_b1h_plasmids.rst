********************
Combine B1H plasmids
********************
[Meng2005]_ uses separate plasmids for the protein "bait" and DNA target 
"prey".  However, this approach doesn't work if both components are libraries, 
because there;s no way to tell which protein went with which DNA target.  For 
my library-vs-library assay, then, I will need to have the protein and the 
target on the same plasmid.

One way to do this is to join the two plasmids in vivo, by expressing a 
recombinase and including recombination sites such that the two barcodes end up 
next to each other in the recombined product.  Gleb and Pierce presented this 
idea in lab meeting about a year ago, so I should talk to them if I want to 
seriously pursue this.  Frankly, though, it seemed complicated and noisy.

The other way to do this is to clone everything into a single plasmid.  This is 
the approach that I want to start exploring here.

Considerations
==============

Asking for help
---------------
It might be worth sending an email to the authors about this, just to ask is 
they have any words of caution.

Plasmid size
------------
Larger plasmids are a little harder to work with and are not transformed as 
efficiently, but I don't think this will be a real problem.  For one thing, I'm 
sure its more efficient to transform one bigger plasmid than two smaller ones.  
For another, the combined plasmid won't be that big.  Assuming I mostly get rid 
of unannotated sequences, it could be as small as ≈6 kb:

  - >2 kb: HIS + URA
  - >1 kb: rpoA + Zif
  - >1 kb: AmpR
  - >1 kb: pSC101 ORI

This would be bigger than most of my plasmids, but not really outside the realm 
of the ordinary.

Expression levels
-----------------
Switching from 2 plasmids to 1 will affect expression levels of the B1H genes, 
because the two plasmids have different copy numbers:

- HIS/URA: pSC101, low copy
- rpoA fusion: p15A, medium copy

I think it's important that the HIS/URA genes be on a low-copy plasmid, so I'm 
going to keep the pSC101 origin and "weak" lac HIS/URA promoter.  This may mean 
that I'll need to increase the expression level of the rpoA-fusion.  Hopefully 
that would just mean adding more IPTG, but I'll cross that bridge when I come 
to it.

Gene insulation
---------------
In the 2-plasmid system, the rpoA-fusion and HIS/URA genes are well insulated 
from expression for any reason other than the intended ones (IPTG induction for 
the rpoA-fusion, rpoA-fusion binding to its target for HIS/URA):

- All of the other genes on the plasmid point in the opposite direction as the 
  rpoA-fusion and HIS/URA, so run-on transcription of those genes will not make 
  transcripts of the B1H genes.

- The HIS/URA gene is preceded by the rrnB T1/T2 terminators.  This would stop 
  transcription from any cryptic promoters pointing in the same direction as 
  the HIS/URA gene, although I don't know if this is really a concern.  These 
  terminators may also serve to prevent rpoA from inducing expression of any 
  genes other than HIS/URA, namely KanR.  I don't know if this is really a 
  concern, either.
  
In making a 1-plasmid system, I'm taking the following steps to maintain as 
much of this insulation as possible:

- Orient the URA/HIS gene in the opposite direction as all of the other genes 
  on the plasmid.  This means orienting the rpoA-fusion in the same direction 
  as other genes.

  I made this decision because I think it's more important to insulate the 
  HIS/URA gene than the rpoA gene.  Expression of the HIS/URA gene *is* the 
  assay, and that expression level may be quite low.  Any increase in basal 
  expression will reduce the sensitivity of the assay.

  In contrast, the rpoA-fusion is controlled by the lac promoter, which is 
  known to be leaky.  So an increase in expression due to run-on transcription 
  is more likely to be negligible.

- Put strong terminators between each gene.  I'll use the terminators from the 
  original plasmids if they're strong, unique, and not preceded by a lot of 
  unannotated sequence.  Otherwise, I'll use terminators from [Chen2013]_ that 
  are strong and have less than 25 bp of identity (to reduce homologous 
  recombination).
  
- Keep the ≈350 bp containing the rrnB terminators between the target and the 
  resistance gene.  Except for the terminators, this sequence is unannotated, 
  but I think it's probably good to keep some space between the target and the 
  nearest non-HIS/URA promoter.

- Keep the rpoA-fusion and HIS/URA promoters far apart from each other.  Even 
  with the ≈350 bp/rrnB insulator mentioned above, I really don't want a 
  feedback loop where rpoA expression induces its own expression.

  This also relates to concerns about getting the barcodes near each other, 
  which is discussed elsewhere.

Reading the barcodes
--------------------
I should think a little about how I'm going to read the barcodes to determine 
which proteins bound which targets.  I have some vague ideas right now, but I 
should talk with someone more experienced to make sure I'm on the right track.

One option may be to not have the barcodes near each other, but to use 
paired-end sequencing to read both of them.  This may not be an efficient use 
of sequencing throughput, but I don't think that'll be limiting anyways.  I'm 
not sure if you can do paired end sequencing with really long fragments (e.g. 3 
kb), though.

Another option is to include restriction sites that can be used to cut out 
everything between the two barcodes.  This could either be a regular enzyme 
(Type IIP) or a Golden Gate enzyme (Type IIS), although the former would 
require a gel purification.

Ultimately, it's premature to think about how I'll do NGS at this point.  I 
will include PCR primers that will allow me to amplify the barcodes, though, 
because I know I'll need that.  And having more PCR primers never hurts.

Target sequence MCS
-------------------
In pH3U3, the target is located in an MCS.  I'm going to keep that as it is 
initially, to simplify the assembly and focus on getting the plasmid to work.  
But ultimately I'll want to replace this sequence:

- An MCS is not really convenient for cloning.  PCR primers and Golden Gate 
  sites are more useful.  In particular, I'll need a PCR primer site on one 
  side of the target to amplify for sequencing.

- Several of the MCS restriction sites are no longer unique in the combined 
  plasmid, so I'd have to go to some effort to make them usable again.

- The MCS is very GC rich.  On one hand this is a problem, because it makes it 
  hard to do PCR, which would be convenient e.g. for changing the target by 
  inverse PCR.  On the other hand, it's conceivable that having a GC-rich 
  region near where rpoA is binding is somehow relevant to B1H.  Something to 
  keep in mind if I have problems.

Restriction sites
-----------------
I'll eventually need to remove all the BsaI and BbsI sites from the plasmid.  I 
can experiment with the plasmid before then, though.

Resistance gene
---------------
I don't like Kan/Chlor, I'm going to switch to AmpR.

Unnecessary elements
--------------------
- I'm going to remove the f1 origin that's present in pH3U3.  As far as I can 
  tell, the f1 origin only makes it so that the plasmid can be propagated in 
  phage (i.e. a "shuttle vector"), and does not play any role in the B1H assay.
Modular cloning
---------------
Since I'm designing this plasmid more or less from scratch, I thought that it 
might be smart to adhere to a modular cloning standard, e.g. MoClo or 
GoldenBraid.  

I just brushed up on MoClo.  The original MoClo system [Weber2011]_ is designed 
for eukaryotic genes, but two E. coli part libraries have been described and 
made available on AddGene.  The first is CIDAR [Iverson2016]_ and the second is 
EcoFlex [Moore2016]_.

- CIDAR mostly use the same overhangs as the original MoClo, but not with the 
  same meanings.  (So MoClo and CIDAR parts are not compatible, but that's 
  fine, they're meant for different organisms anyways.)  It's not clear to me 
  how CIDAR transcriptional units (TUs) are assembled, but presumably I'm just 
  missing something.  AddGene has both a CIDAR kit and a CIDAR extension kit, 
  which total to more parts than EcoFlex has.

- EcoFlex uses completely different overhangs than MoClo.  It also has support 
  for N-terminal tags.

- I don't think I can directly use CIDAR/EcoFlex, because I want my genes 
  pointing in opposite directions.

Ultimately, I think my goals are different enough that using a modular cloning 
system would be more effort than it would be worth.  These systems are really 
meant to facilitate pathway engineering.  In this context, it's important to be 
able to create assemblies with many genes, and to be able to easily try 
different promoters for each gene.  I don't have that many genes or that many 
promoters/parts to try.  I also have concerns like getting good transformation 
efficiency or being able to read my barcodes that these modular systems may get 
in the way of.

Assembly
========
6 parts:

- HIS/URA from pH3U3
- rpoA-zif268 from pB1H1
- AmpR from p004
- pSC101 from pH3U3
- Terminators from gBlock
