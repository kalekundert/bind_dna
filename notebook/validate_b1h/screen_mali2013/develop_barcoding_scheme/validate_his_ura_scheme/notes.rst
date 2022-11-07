***********************
Validate HIS/URA scheme
***********************

This scheme is a modified version of the "single-barcode scheme" described in 
:expt:`180`.  The only difference is that, since we're using HIS/URA as a 
reporter instead of RNA seq, there's no need to express the barcode.  (In fact 
it's better not to.)  This ends up changing how the cloning would be done.

Considerations
==============
- Because I want the whole DBD oligo to be a CDS, and I don't want to express 
  the barcode, I'm going to need to insert either a promoter or terminator 
  sequence between the CDS and the barcode.  There's not gonna be any way to 
  assemble the library without extra cloning steps.

- If I keep the barcode upstream of the target, I could maybe get away with 
  including a terminator in the target oligo and skipping a cloning step.  But 
  it's probably cheaper to make the oligos as short as possible, and to install 
  as much constant sequence as feasible with regular cloning steps.

- If I'm going to clone constant sequences into both the target and the DBD 
  anyways, it doesn't really matter where I put the barcodes.  So I should just 
  choose a site that I would least expect to interfere with anything.  My gut 
  says after the His/Ura cassette.

  - I could measure this by creating libraries with random barcodes, but only 
    one target/DBD and seeing which gives greater variability.

- I could keep the barcodes separate, then use a cloning/recombination step to 
  bring them together after the selection.  But I think I'd rather do a little 
  bit of extra cloning in advance than complicate the sequencing step at the 
  end.

  - There's actually an established technique for doing this, called `mate-pair 
    sequencing`__.  Basically, the process goes: (i) amplify a 
    several-kilobase-long fragment, (ii) ligate tags onto the ends, (iii) 
    circularize, (iv) re-fragment, (v) ligate adaptors, and (vi) sequence.  
    Those reads with the tags in the middle contain both ends of the initial 
    amplicon.  But it's still possible for intermolecular ligation to occur, 
    and some fraction of the read depth is wasted on constant sequences from 
    the interior of the amplicon.  Overall, I still think it's better to do the 
    cloning in advance and have something simple to work with at the end.

__ qute://pdfjs/web/viewer.html?filename=tmp4sq5d9e__technote_nextera_matepair_data_processing.pdf&file=&source=https://www.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

- I never did the experiments to see which parts of the plasmid can be 
  truncated.  I planned such experiments for the p194 scaffold, which wouldn't 
  necessarily apply to this scaffold anyways.  I did try to add back the 
  terminator after the rpoZ-DBD (:expt:`161`), but that wasn't successful.

  For now, I should move forward with the p189 scaffold.  But I might try to 
  restart the truncation experiments while I'm in the process of cloning the 
  test library.

- Stefan, Tina, and Erik tell me that restriction cloning is the most efficient 
  way to assemble a library, so I should use Golden Gate assembly instead of 
  Gibson assembly.  This also helps avoid PCR steps, which can introduce errors 
  and biases.

- I'll have to remove a number of existing BsmBI and BbsI sites from my 
  sequences.  I might just order gBlocks to get around the cloning steps.

Oligo-pool PCR
--------------
- Twist protocol: https://tinyurl.com/3sce5397
- KAPA Hi-Fi protocol: https://tinyurl.com/2p8xdzdm

- You want just enough cycles to get the amount of material you need.  In each 
  round of amplification, different templates may be amplified unevenly based 
  on GC-content, particular sequence motifs, etc.  The fewer cycles there are, 
  the less these biases will be magnified.

  The KAPA manual claims that 100 ng - 1 µg is usually enough, although of 
  course it depends on the application.

- Tina recommends starting with a small number of cycles, running 1 µL on an 
  E-gel to see if the amplicon is just barely visible, then adding cycles until 
  it is.

  The detection limit of E-gel EX is ≈1 ng, so this implies a final 
  concentration of 1 ng/µL and a final yield of 25 ng (for a 25 µL reaction).  
  This is pretty much in line with the KAPA recommendations.

- KAPA primers should include a 3' phosphorothioate modification, to prevent 
  the strong proof-reading activity of the KAPA polymerase from degrading the 
  primers.

- Over-amplification: the KAPA protocol has a good description:

  In library amplification reactions (set up according to the recommended 
  protocol), primers are typically depleted before dNTPs. When DNA synthesis 
  can no longer take place due to primer depletion, subsequent rounds of DNA 
  denaturation and annealing result in the separation of complementary DNA 
  strands, followed by imperfect annealing to non-complementary partners. This 
  presumably results in the formation of so-called "daisy-chains" or tangled 
  knots, comprising large assemblies of improperly annealed, partially 
  double-stranded, heteroduplex DNA.
  
  These species migrate slower and are observed as secondary, higher molecular 
  weight peaks during the electrophoretic analysis of amplified libraries.  
  However, they are typically comprised of library molecules of the desired 
  length, which are separated during denaturation prior to target enrichment 
  (using hybridization capture protocols) or cluster amplification.  Since 
  these heteroduplexes contain significant portions of single-stranded DNA, 
  over-amplification leads to the under-quantification of library molecules 
  with assays employing dsDNA-binding dyes.

Purifying library inserts
-------------------------
Initially I was planning to PAGE-purify the library oligos after amplifying and 
digesting them.  However, doing so would lose a lot of material.  I would need 
to do more PCR to make up for that, which in turn would introduce more errors 
and more bias.  It might not be worth it.

The purpose of the purification is to get rid to the sticky ends that I don't 
want to participate in the ligation reaction.  For the constant regions, a gel 
purification is necessary, because the plasmid will be cut into two parts and I 
only want one.  But I can make as much of those plasmids as I want, so it 
doesn't matter if I lose some.

My options for the library inserts:

- No purification: Hope that the ends produced by the digest don't 
  significantly interfere with the ligation reaction

- Column purification: Not viable for l3, because the product of interest is 
  only 59 bp (compared to 26 bp for the cleaved-off ends).

- Bead purification: Might allow elution in smaller volumes, and can (in 
  theory) tune bead ratio to get rid of digested ends.  But again, probably not 
  viable for l3.

- PAGE purification: Will likely have to do more PCR cycles.

- Streptavidin pull-down: I could amplify with biotin-modified primers, then 
  use streptavidin beads to remove the digested ends.  This would require fancy 
  reagents (although I already have the beads), and wouldn't remove ends that 
  weren't synthesized by PCR, but could give good purification with little 
  loss.

I should probably try these and see which works the best, where "working the 
best" would be defined as:

- Most transformants
- Most correct sequences (by NGS)

The number of transformants won't really matter unless it's awful, because the 
sub-libraries are small and I should be able to get full coverage even if the 
cloning isn't perfectly optimized.  The real interesting part will be to see 
which gives the best fidelity.

To begin with, though, the simplest thing is to just not do any purification at 
all, so I'll start with that.

Approach
========
.. figure:: his_ura_cloning.svg

Reaction #0: pUC backbone
-------------------------
- Cleavage near ends:

  - HindIII: 3 bp
  - NotI: >5 bp

Reaction #1: DBD + pSC101
-------------------------
- This needs to be scarless, because the constant region could include the 
  C-terminal end of the DBD.  That means Gibson or Golden Gate, and as 
  described above, Golden Gate is preferred for library assembly.

- Cleavage near ends:

  - BsaI: 1 bp
  - BsmBI: 1 bp

Reaction #2: DBD + rpoZ + pSC101
--------------------------------
- This needs to be scarless because I'm making a fusion protein, and I might 
  want to include more or less of the DBD in the constant region.

- Put constant sequence in a plasmid, so that I don't need to do PCR to add 
  BsmBI sites.  The primers would have to be very long, and so the chance of 
  errors would be high (not even counting errors due to the polymerase).  
  Plasmid amplified DNA can be sequenced and propagated with high fidelity.

- I should be able to mix the backbone and the library fragment in one pot with 
  the restriction enzyme, but I also suspect that I might get better 
  transformation efficiency if I digest and purify both fragments separately.  
  It's probably worth doing an experiment to see what works better.

Reaction #3: target + pUC
-------------------------
- The sequence upstream of the target very GC-rich and not amenable to PCR.  I 
  currently use restriction cloning when working with this sequence.

- I need to make a pUC vector with a terminator and a NotI site

- Restriction sites:
  - NotI (I'd need to remove the NotI site from the FLAG sequence)
  - Something from the pUC MCS, e.g. BamHI
  
- I'm going to gel purify the backbone after the restriction digest, so it 
  might be helpful if I can visually distinguish single-cut and double-cut 
  plasmid.  This means I want an insert that is as big as possible, and a 
  backbone that is as small as possible.

  - The backbone will be about 2.4 kb.  I could shrink that a little, but not 
    below 2 kb.

  - The insert could be any length.  I guess I'm only limited by the fact that 
    I want to buy it as a gBlock, so I don't want to spend too much.  The 
    gBlock in question will have ≈400 bp devoted to functional sequence.  The 
    rest can just increase the size of the insert.  All gBlocks up to 500 bp 
    cost the same amount, and after that the price scales linearly.

    I'm going to start with a 500 bp gBlock, meaning a 100 bp insert.  That  
    might not be easy to see, but it's cheap and will likely work well.  If 
    necessary, I can make the insert longer with a single cloning step.

- Cleavage near ends:

  - BsaI: 1 bp
  - BsmBI: 1 bp
  - NotI: >5bp
  - BsrGI: 1 bp
  - HindIII: 3 bp

- Sequences:

  - gtgtacacccgggcggccgctGCGTGGGCG 

    - BsrGI, NotI, target site (uppercase)
    - Include all of BsrGI because I don't know how much distance from end is 
      necessary to get full cleavage, and I'd rather be safe.   I'm also not 
      going to be limited by the length of this oligo.

  - ggactGAGACG

    - Seamless junction (ggac), BsmBI site (uppercase)
    - Use T as spacer because current sequence is very GC rich, and there are 
      A's nearby but not T's.

Reaction #4: target + His/Ura + pUC
-----------------------------------
- I've avoided doing PCR near the target site, but that's because of the really 
  high GC content upstream of the target.  Downstream of the target is the lac 
  promoter, and while it's a bit AT-rich, it has plenty of viable primer 
  sequences.  So I think Gibson is the way to go here.

  - pUC has the same lac promoter as my gene.  That's slightly annoying.  It 
    means I'll want to cut the lac promoter out of the pUC plasmid when adding 
    the terminator etc.

- Need to make plasmid for constant fragment:

  - Avoid using PCR to generate fragment.  PCR has higher error rates than 
    miniprep, and in this case would depend on the fidelity of long primers.

  - Need to cure two BsmBI sites.

    - Delete the one in the f1 ori.
    - Silently mutate the one in URA3.

  - Cloning strategies:
    
    3-step:

    - Gibson to put HIS/URA cassette in pUC backbone
    - PCR to delete site in f1 ori (delete a lot, or just a little?)
    - PCR to silently mutate site in URA3

    1-step:

    - 3-part Gibson assembly with:

      - pUC backbone
      - HIS/URA amplified from existing plasmid
      - URA/terminator ordered as ≈1kb gBlock ($230)

    Inclined to save time and do the assembly.

Reaction #5: final assembly
---------------------------
