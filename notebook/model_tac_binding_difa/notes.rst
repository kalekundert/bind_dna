**********************
Model Tac binding difA
**********************

Ben is interested in having a plausible model of Tac binding its target site 
(difA), based on the crystal structure of Cre binding loxP.  Tac and Cre are 
homologous.  Both the N- and C-terminal domains superimpose reasonably well, 
but that the linker is in a different conformation.  Below is the specific 
request from Ben:

   This is something that I'm thinking one thing that we can do first is to try 
   to properly dock this crystal of Tac (PDB: 5HXY) onto its target site difA 
   (TCGCTGTATAActatcgTTATACGGCAT) using Cre pre-cleavage synapse complex as a 
   guide (PDB: 1Q3U … or 2HOF?). I have aligned the C-terminal domains with 
   each other pretty well, but I have not delved into aligning the N-terminal 
   domains.  It will be good to know which residues are base-contacting. I can 
   sort of guestimate with the alignment I have now, but it would be good to 
   have a more solid picture. Also, would be good to test the chimeric proteins 
   on this model as well.

I actually don't need to do any docking, because it's clear that both the N- 
and C-terminal domains should be superimposed on each other.  As I see it, 
there are two things that modeling may be useful for:

- Determine plausible atomic interactions between Tac and difA.  This will 
  basically mean installing the mutated sequence and relaxing with backbone 
  constraints.

- Determine the conformation of the linker connecting the two DNA-binding 
  domains.  In particular, one of the helices in the N-terminal domain is 
  unraveled in the Tac structure; presumably this helix should be folded.

After talking with Ben, I decided to go ahead and model the atomic interactions 
between Tac and difA, i.e. relax the superimposed domains.  This could 
influence decisions about which mutants to test in the recombination assay.  
The linker conformation is less important.

Notes
=====
- If I made mistakes in the DNA, that can mess up the backbone.

- Repacking is not enough to get reasonable base pair conformations, but 
  relaxing seems to work well.

Protocol
========
1. Split :pdb:`5HXY` into two domains:

      > load 5hxy.pdb
      > select nterm, resi 4-97
      > save 001_n.pdb, nterm
      > save 001_c.pdb, not nterm

2. Superimpose both domains onto chains A and B of 1q3u (Cre):

      > load 1q3u.pdb
      > load 001_n.pdb
      > load 001_n.pdb
      > split_states 001_n
      > load 001_c.pdb
      > load 001_c.pdb
      > split_states 001_c
      > super 001_n_0001, 1q3u and chain A
      > super 001_c_0001, 1q3u and chain A
      > super 001_n_0002, 1q3u and chain B
      > super 001_c_0002, 1q3u and chain B
      > remove resn HOH+MG+IOD
      > save 002_super_on_cre.pse
      > remove 1q3u and not chain C+D
      > save 003_domains.pdb

   - I had the thought that maybe I should've removed any residues that 
     obviously differ between the two crystal structures (e.g. the linker loop) 
     before doing the superposition, so that the algorithm wouldn't try to 
     superimpose them.  However, the ``super`` algorithm does this 
     automatically, so there's no problem with superimposing the domains like I 
     did.

   - The ``align`` algorithm in pymol creates a sequence alignment and performs 
     a superposition based on that alignment.  The docs claim that it works 
     best with sequence identity >30%.  In contrast, the ``super`` algorithm in 
     pymol is purely geometry based (it doesn't consider sequence at all) and 
     is better for structure with low sequence similarity but high structural 
     similarity.  To decide which algorithm to use, I performed a sequence 
     alignment of Cre and XerA using the `EMBL Needleman-Wunsh web server`__.  
     That gave a sequence identity of 17.6%, suggesting that the ``super`` 
     algorithm is more appropriate.
     
     __ https://www.ebi.ac.uk/Tools/psa/emboss_needle/

   - I removed all the water molecules.  Practically, this is because the 
     waters didn't superimpose with the rest of the protein, and it seemed like 
     it'd be a pain to get them to.  I looked at some of the waters with low 
     B-factors, but they didn't seem particularly important.  And of course, 
     none of the interfacial waters would be in this structure.  That will be a 
     caveat of the model, though.  Any waters in the interface will not be 
     predicted.  That could also lead to incorrect rotamer predictions.
   
   - I also removed the Mg and I atoms associated with the DNA, because they 
     didn't seem to be forming specific interactions.

   - Note that ``003_domains.pdb`` only appears to have one monomer when opened 
     in pymol.  This is because the two monomers have exactly the same atom, 
     chain, and residues numbers.  Show lines and you'll see everything.

3. Manually edit the PDB file to put the domains in the correct order and to 
   give the protein chains different names (A and B).  I decided to keep the 
   native residue numbers, because I think it'll make the model easier for Ben 
   to use in the end.  I'll just deal with any extra complexity this causes for 
   Rosetta.

4. Manually edit the PDB file to change the loxP sequence to the difA sequence.  
   For nucleotides that are the same in both sequences, I didn't do anything.  
   For nucleotides that differed between the two sequences, I deleted the 
   non-backbone atoms and changed the residue name of the remaining atoms.

   - Change UMP nucleotide (uridine) to DT (thymidine).  Not sure why this was 
     there in the first place, although the crystal structure paper probably 
     explains.

   - Chain C is the top strand::

               100       110       120       130
               0123456789012345678901234567890123456789
      difA:         TCGCTGTATAActatcgTTATACGGCAT
      loxP:      ATAACTTCGTATAatgtatgcTATACGAAGTTAT
      Chain C: CGATAACTTCGTATAATGTATGCTATACGAAGTTATC

               100       110       120       130
               0123456789012345678901234567890123456789
      difA rc:      ATGCCGTATAAcgatagTTATACAGCGA
      loxP rc:   ATAACTTCGTATAgcatacatTATACGAAGTTAT
      Chain D: GGATAACTTCGTATAGCATACATTATACGAAGTTATC

   - Make sure I didn't make any mistakes::

      $ ./check_resn_per_resi.py 005_fix_difA.pdb

5. Manually prepare an input model for relaxation:

   - If I'd been thinking ahead, I could've done all this is the first few 
     steps.  But that's fine.

   - Delete the linker between the N- and C-terminal domains (residues 87--108) 
     from the model.  I'm going to run ``dbp_relax_b``, which is going to 
     optimize the restraints such that 90% of B-factors are satisfied, so I 
     want to exclude these atoms because I know they aren't in the right place 
     (so they should be allowed to move more.

   - Delete the C-terminus of chain B (residues 266--283) from the model.  This 
     is part of the dimer interface, and it's clearly in the wrong 
     conformation.  I'll probably want to do a 

   - Put the two domains in separate chains.  I don't expect the chains to move 
     much, but I want them to move freely and without regard for any weird 
     fake-bonds from the deleted residues connecting the domains.
     
6. Relax the model::

      $ dbp_relax_b init 007_relax_XerA 006_unrelaxed.pdb

   Edit ``007_relax_XerA/conf.toml``.  Change ``preoptimize.percent-within`` 
   and ``optimize.percent-within`` to 75, to account for the fact that this 
   model will have more atoms out-of-place than a model built directly from a 
   crystal structure.

7. Connect the N- and C- terminal domains via loop modeling.

   .. update:: 2019/07/31

      I didn't actually do this step.

   - I've gotten curious about this linker because it doesn't really seem long 
     enough to connect the two domains.  I mean, it fits, but it seems really 
     stretched out.  This linker also seems to replace a whole helix in Cre.
     
   - fKIC would probably be the best way to make a prediction about what that 
     loop would look like:
     
      - I don't know whether or not there should be secondary structure in the 
        linker, but if there should, fKIC would be more likely than anything 
        else to find it.
        
      - This is also a long loop (~20 res), and fKIC performs much better than 
        CCD or NGK on long loops.

   - Download a fasta file for XerA from `Uniprot`__. :download:`xera.fasta`

     __ https://www.uniprot.org/uniprot/Q9HIM5

   - Submit the above fasta file to the `Robetta fragment server`__.

     __ http://robetta.bakerlab.org/fragmentsubmit.jsp

     Note that it's important to use the fasta file from UniProt rather than 
     one generated from ``5hxy.pdb``, because the PDB file is missing some 
     residues on the N-terminus.  These residue should be considered when 
     making fragments, and leaving them out causes indexing problems.

Results
=======
- I allowed very little backbone movement.  Basically I assumed that the 
  :pdb:`5HXY` crystal structure and the alignment onto :pdb:`1Q3U` were more 
  accurate than anything Rosetta would come up with, so I didn't allow Rosetta 
  to deviate much from that template.  Of course, this means that the model is 
  incapable of predicting any major conformational changes.

- This model does not include explicit water molecules.  This is definitely a 
  weakness, because water-mediated DNA interactions will not be accurately 
  predicted.

- I didn't model water-mediated H-bonds at all.  These H-bonds are quite common 
  in DNA interfaces, so some of the direct DNA contacts in the model may be 
  water-mediated in real life.

- I gave the backbone very little freedom to move.  Basically I took the 
  backbone structure from 5HXY and superimposed the N- and C-terminal domains 
  onto 1Q3U.  I didn't trust Rosetta to be more accurate than that, so I didn't 
  let it move the backbone much more than ~0.5A.  It's possible that the real 
  XerA is differs from Cre in ways this model wasn't allowed to consider.

- The model doesn't include the linker between the N- and C-terminal domains, 
  or one of the helices involved in the dimer-dimer interface.  The 
  conformations that these regions adopted in the 5HXY structure were clearly 
  incompatible with the DNA-bound model, so I left them out.  It would be 
  possible to model these regions, too, but I didn't think the results would be 
  interesting enough to merit the effort.

- Adenine 118 is flipped out of the double-helix in the model.  This is 
  obviously wrong.  The cause is probably related to the fact that the DNA 
  backbone in really curved in that region, and maybe loxP and difA curve 
  slightly differently.  In any case, I don't think this specific error affects 
  the any of the protein-DNA interfaces, so I'm not too worried about it.

- I checked in pymol to see (qualitatively) how well the electrostatic surface 
  of the model corresponded to the DNA interface.  The correspondence was good, 
  but not as good as 1q3u:

   - E28 and E210 both seem to be making relatively unfavorable contacts with 
     the DNA.

   - I suspect that the linkers between the N- and C-terminal domains help bind 
     the backbone, but the model doesn't include those linkers.

- Final disclaimer: This is a pretty big model, and I don't really know how 
  well Rosetta was able to handle it.  I tested the algorithm on a smaller 
  protein (a Zn-finger) and it qualitatively seemed to give reasonable results.  
  But it's possible that the algorithm wasn't able to converge as well on such 
  a big system.  To really know how well the algorithm could be expected to 
  work, I'd have to benchmark a bunch of structures of similar size for 
  Cre/XerA, which would be way too much work for something like this.  So just 
  take everything with a grain of salt.
      
References
==========
Protein-protein docking tutorial:

https://www.rosettacommons.org/demos/latest/tutorials/Protein-Protein-Docking/Protein-Protein-Docking
