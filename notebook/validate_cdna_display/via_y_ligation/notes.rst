**************
Via Y-ligation
**************
The cDNA-display protocol from [Naimudden2016]_ uses y-ligation, and claims 
that it gives higher efficiency than the more traditional splint-ligation.  
Here I'm going to try the [Naimudden2016]_ protocol, with some 
PURExpress-specific details taken from [Barendt2013]_.

.. toctree::
   :hidden:

   order_linker_n/notes
   make_linker_via_click_chemistry/notes
   linearize_cdna_display_gene/notes
   optimize_mrna_transcription/notes
   eliminate_mrna_degradation/notes
   eliminate_mrna_linker_degradation/notes
   validate_pseudo_linker/notes
   denature_unligated_linker_n/notes
   anneal_linker_to_mrna/notes
   ligate_linker_to_mrna/notes
   remove_unligated_linker/notes
   remove_unligated_mrna_via_page/notes
   remove_unligated_mrna_via_biotin/notes
   express_mrna_via_purexpress/notes
   express_mrna_via_purefrex/notes
   express_mrna_via_nebexpress/notes
   express_mrna_via_s30_lysate/notes
   attach_mrna_to_mwasabi/notes
   attach_mrna_to_flag/notes
   attach_mrna_to_zif268/notes

Considerations
==============

C-terminus
----------
[Naimudden2016]_ calls for the following sequence elements C-terminal of the 
gene of interest:

- 2x GGGS
- His-tag
- 1x GGGS
- Y-tag

I'm curious if these C-terminal elements can't be improved.

- Stop codon with ΔRF PURExpress

   - The idea with using a stop codon is that, if there are no release factors 
     in the reaction, the ribosome will stall on the stop codon, which will the 
     puromycin a chance to react.

   - I might even be able to subsequently add release factors to help release 
     the puromycin-linked protein from the ribosome.  None of the papers I've 
     read have seemed to say much about releasing from the ribosome, so maybe I 
     don't need to worry about it.  But if it turns out to be a problem, this 
     could be a solution.

   - Using a stop codon would require me to use PURExpress, because any lysate 
     would include release factors.

   - A stop codon would make the C-terminus more consistent.  Currently I 
     believe that translation will stop when the ribosome encounters the 
     double-stranded region where Linker-N is bound, but it's also possible 
     that Linker-N itself will be translated (it's mostly Gly).  A stop codon 
     would make it more clear where translation stops.  
     
     Unfortunately, though, I don't really have a way of knowing how consistent 
     the C-terminus is.  So I won't know if this is really a problem, or if a 
     stop codon even does anything to help.  If I really wanted to know this, 
     I'd probably have to do mass spec.

   - A stop codon would also enable the use of 3' barcodes.  See :expt:`53`.

   - The poly-A arm supporting the puromycin in Linker-N is optimized for use 
     without a stop codon.  Granted, it's also optimized for use with 
     eukaryotic ribosomes, which are slightly bigger than bacterial ones.  
     Still, adding a stop codon would move the puromycin further from the 
     A-site, which may interfere with the ability of puromycin to attack the 
     stalled peptide.

  All told, using a stop codon could be a good thing or a bad thing.  It's not 
  something that [Naimudden2016]_ would've tried, because they were using 
  rabbit reticulocyte lysate.  So I think it's worth trying, but I should get 
  the system working first without it.

- Remove or shorten the GGGS linkers

   - First of all, [Naimudden2016]_ used GGGS linkers, while most everyone else 
     uses GGGGS (4 Gly instead of 3).  It probably doesn't matter, but I'm 
     tempted to use the more canonical (longer) linkers.

   - I'm skeptical that 2 repeats of the Gly/Ser linker are needed before the 
     His-tag.  After all, the His-tag and the second linker can be thought of 
     as one big linker between the protein and the mRNA/cDNA.  Again, though, 
     this probably doesn't matter too much.

Variable transcript length
--------------------------
As described very clearly by [Gholamalipour2018]_, T7 polymerase can add a 
variable number of extra nucleotides to the end of transcripts.  This is caused 
by the transcripts folding on themselves and priming additional transcription.   
This phenomenon occurs more readily when the concentration of RNA is high (as 
in IVT reactions), and in these conditions [Gholamalipour2018]_ show that the 
correct-length transcripts can become a minor product.

It's possible that this extra transcription could interfere with the ligating 
of Linker-N to the 3' end of the transcript.  However, the actual Y-tag 
sequence seems like it would be fairly resilient to this effect, because it it 
mostly G and therefore has poor self-complementarity.  Also, [Naimudden2016]_ 
reported good ligation yields, and they didn't mention this issue at all.  So 
I'm not going to worry about this for now, but I wanted to at least put my 
thoughts down.
