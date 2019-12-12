******************
Ideas and Feedback
******************

2019/04/03:

From talking with Gleb:

- Think about ways to fail quickly.

   - I think I'm pretty good on this.  It should be pretty easy to progress 
     through testing:

      - One protein vs one target.
      - One protein vs one target in pool of decoys.
      - One protein vs targets with different affinities.
      - Repeat all above tests with different protein classes.
      - Multiple proteins vs multiple targets.

- Make sure I'm not blocked if I can't get the assay to work.

   - I don't really doubt that the assay will "work".  My bigger concern is 
     that the throughput and signal-to-noise ratio might not be good enough to 
     do truly high-throughput screening.

   - Some things I should keep in mind if I'm having trouble:

      - There are people in the lab using droplet emulsions for library 
        screening.  Maybe I could do something like that.

   - In the case that I can't get the cDNA-ligation assay to work, I can:

      - Try other assays, like [Gu2014]_.

      - Revert to testing small numbers of designs using specificity assays 
        like [Zykovich2009]_ or maybe [Meng2007]_.

- Test a variety of control scaffolds.

   - One, it is important for different scaffolds to work if I want to talk 
     about general-purpose design algorithms.

   - Also, it's possible that some scaffolds will just work and others just 
     won't.  For example, the orientation between the C-terminus and the 
     binding interface could affect the efficiency of the assay.  This is 
     something I want to get a feel for.

2019/04/08:

From talking with David Johnson:

- There are definitely reasons to prefer display technologies:

   - Optimization is much easier; don't need to mess with promoters, RBS, etc.

   - Throughput is much higher; not limited by transformation efficiency.

   - Control of conditions is better; pH, salt, temperature, etc.

- Consider an emulsion method:

   - Protocol: 

      - Link barcoded protein genes to target DNAs enzymatically.
      - Encapsulate linked DNA in emulsion droplets.
      - Optionally cut the gene off of the target (leaving the barcode).
      - Express protein.
      - Crosslink the protein to the DNA.
      - Break the emulsion.
      - Pull down the proteins (e.g. via His-tag).
      - Uncrosslink
      - Sequence

   - Not sure how hard it would be to add things (e.g. IVTT mix, restriction 
     enzyme) to droplets like this.  Would it require custom microfluidics?

   - Crosslinking might make it hard to distinguish weak vs strong 
     interactions.

   - Formaldehyde crosslinks are reversible.  See the CHiP-Seq literature for 
     details.

- Consider low-throughput methods:

   - This is how big pharma does their screens: no fancy method development, 
     just brute force.  Obviously we can't brute force as much as they can, but 
     we might be able to brute-force enough to answer our question.

   - Potential advantages:
      
      - Less messiness in terms of cross reactivity, dimerization, etc.

      - Less method development work required.

      - Should still be feasible to test ~1000 designs, which may be be enough 

   - A good assay might be based on [Zykovich2009]_ (and actually, this is 
     similar to the emulsion protocol described above):

      - Order designs:

         - This would be expensive, because I'd have to order individual 
           gBlocks, I think.  A full 96-well plate would be >$10,000, which is 
           probably affordable.  I might not be able to get to 1000, though.

         - Less limited in terms of size (compared to oligo pool).  I'd still 
           be limited by IVTT, though.

         - Maybe I could order oligo pool, transform, and pick colonies.  I'd 
           have the problem where I'd have to test ~10k designs to be sure I 
           got my whole library of 1000, and I might have to do ~10k minipreps, 
           but it might be more plausible price-wise.  I wonder if I could just 
           lyse the cells and add IVTT mix, rather than purifying the 
           plasmids...  But wait, I'd need some kind of sequencing step to know 
           which protein went with which barcode...

         - I might also be able to do an oligo pool with a bunch of subpools 
           (e.g. amplify with orthogonal PCR primers).

      - Clone protein library with a pulldown tag.
         
         - [Zykovich2009]_ used MBP/amylose.
            - Pros: probably helps with expression and stability.
            - Cons: MBP is >400 aa, maybe prohibitive in terms of in vitro 
              expression.

         - Don't see why I couldn't use His-6 or FLAG instead.

      - Express protein in vitro in in 384/1536-well plates.

         - Could probably use pipetting robot to improve throughput.

      - Prepare 384/1536-well plates with DNA barcodes that correspond to the 
        proteins.

      - Prepare 384/1536-well plates with target library in each well.

      - Attach the barcodes to the targets:

         - Make compatible sticky ends.
         - Add ligase.
         - Denature everything.

      - Add protein, incubate.
      - Crosslink (formaldehyde)
      - Quench (glycine or tris)
      - Pulldown
      - Reverse crosslink (high salt)
      - Sequence

- Thoughts abouts cDNA-display:

   - The problem with mRNA- and cDNA-display is that it's hard to get it to 
     work if you're not being taught by someone who really knows how to do it.  
     In other words, it's something of an art.

      - I really don't need good yield, though, so I might be able to get it to 
        work well enough pretty quickly.

   - Knowing how to do mRNA- and cDNA-display could be a valuable skill in 
     later stages of my career.  Since I'm at the beginning of my postdoc, if 
     it's something I want to have in my repertoire, it might be worth 
     investing the effort to learn now.

      - I don't really think of myself as someone who's going to want to pound 
        out massive libraries.  But I do think that DNA is a great readout, and 
        so attaching DNA to things could be a generally useful approach.

- Thoughts about my ligation assay:

   - Problem: Protein might bind to its own cDNA, because tagged cDNA has 
     extremely high local concentration.

      - Might end up seeing things that don't bind very tightly to DNA, because 
        the things that do bind tightly will be bound to their own cDNA.

      - I might be able to mitigate this by adding higher levels of target, but 
        the local concentration of cDNA will always be the highest.

   - Could maybe try using stop codons without release factors to get 3' 
     barcode.

      - There is a PURExpress kit with the release factors eliminated.

      - This is something I should probably try.

   - My thought at the end of the discussion is that it's worth giving my 
     ligation assay a try.  Of all the assays I'm considering, it probably has 
     the most upside if it works.  And I should be able to get an idea of 
     whether or not it will work pretty quickly.


From talking with Kettner:

- Thoughts for my assay:

   - Pooled all-by-all library: May take a long time to reach equilibrium, will 
     need to try longer and shorter incubation times to see when results 
     stabilize.

- Consider "cis-display"/"cad-display"
   
   - Protocol:

      - Fuse RepA to C-terminus of protein of interest.

      - RNA polymerase stalls on CIS sequence (terminator?)

      - RepA:

         - Nucleates folding on CIS sequence?
         - Then binds adjacent ORI?
         - Non-covalent, but apparently strong and long lived...
         - 33 kDa

      - Ta da.

   - Advantages:

      - Easy to try, no custom reagents
      - Proof-of-concept?

   - Disadvantages:

      - Noncovalent
      - >300 aa.

- Disadvantages of cDNA-display:

   - Can get partial products.

      - Not too concerned about this.
      - Can purify with C-terminal tag.
      - My proteins are short enough (and my libraries small enough) that this 
        shouldn't really be a problem.

   - Time consuming steps:

      - Overnight ligation.

      - Purifications:

         - Purify for successfully ligated mRNA/DNA:

            - Two purifications

               - One for mRNA
               - Second for DNA linker

            - Don't want free puromycin linker:

               - Will cause aborted translation.

            - Don't want too much unlabeled mRNA:

               - If 90% of protein is unlabeled, assay might not work well
               - Wasting IVTT capacity.
               - May be able to get away without this purification if ligation 
                 is efficient.

         - Purify for fully expressed protein.

      - Not too concerned about it taking a few days to make display libraries, 
        though.

         - I'm not doing multiple rounds of selection; I'm just doing the 
           library prep once for each experiment.
         - I don't think the assay will be the bottleneck.  It will take longer 
           to generate designs, analyze the results, and tweak the algorithms.

- Alternative binding assay:

   - Ligate targets to protein genes.

   - Create DNA-display library (i.e. either cDNA-display of CIS-display)

      - Need to get good display library efficiency, because now I need to 
        express ~1e10 different transcripts.  But that is still within the 
        limits of what in vitro display should be capable of.

   - Allow protein to bind target on it's own DNA.

      - Each protein is only supposed to bind its own target.
      - Keep the reaction dilute to ensure that interactions are predominantly 
        intramolecular.
      - Can add free target to increase stringency for binding.

   - Incubate with DNase.

      - Binding should confer protection from nuclease treatment.
      - Can tune intensity of nuclease treatment.

   - Sequence to identify targets that were protected from DNase treatment.

- To get cDNA linker, just order from a company: Midlands something?  Midlands 
  CRC?  Maybe ask Kettner for details.

2019/04/09:

Talk with Gabe:

- Talk to Rhiju

- Think about how to categorize results, so I can find metrics that correlate 
  with success.

2019/04/29:

- It might be useful to make my cDNA be circular on one end.  This would impart 
  a directionality to ligation and a resistance to exonuclease degradation.  I 
  could do it by ordering a small hairpin as a cap, then using a restriction 
  enzyme and ligase to attach it to my DNA where I want it.

2019/05/15:

Feedback from Jacob Corn:

- Design of recombinases is very exciting. Definitely a high bar, but very 
  worthwhile. George's lab has been trying to do this for a long time, and it's 
  just as important now as BC (before CRISPR). People stopped in part b/c it's 
  so hard, and CRISPR is so easy. There are still good applications for 
  protein-DNA interface design, but many people have decided that it's easier 
  to redesign Cas9 (or jam effectors onto Cas9) than it is to design a 
  protein-DNA interface. They are not necessarily right.

- Recombinase design is much harder than meganuclease design, since MNs don't 
  need to move.

- MN design was pretty easy, but still monstrously hard. Several of the MN 
  successes were accidents based more on brute force than good design. One 
  problem is that protein-DNA interfaces are far more wet than protein-protein 
  interfaces. Water bridges are quite common, and in fact almost as common as 
  direct hydrogen bonds. See classic papers from Janet Thornton's lab (first 
  author Nick Luscombe). These came out in the early 2000s so could be updated.  
  But stress the importance of solvation.

- Water is quite hard, even when using solvated rotamers. There is a reason 
  that MD people have so many water models. Phil Bradley has spent a lot of 
  time thinking about this. Might be worthwhile to pick his brain.

- I like the lib-on-lib approach to test a lot of designs. Definitely key for 
  success. Keep in mind that gaining knowledge from the successful designs is 
  do-able. But it's almost impossible to gain knowledge from unsuccessful 
  designs. To paraphrase Tolstoy, each bad design is bad in its own way. It 
  takes a lot of work to debrief bad designs and learn why they failed and 
  hence learn from them. It's still good to talk about this in a proposal 
  (since people outside the field like to hear it as a motivation). But don't 
  have high hopes, especially with protein-DNA design.

What I distilled from the above feedback:

- Solvation and catalysis:
  
   - I like the idea that for each design task, there's one big *thing* you 
     need to get right:
      
      - Backbone H-bonds for structured loops
      - Glycines to relieve strain in β-barrels

   - Of course, this doesn't have to be true.  Any design will need to get lots 
     of things right.  What I'm assuming is that most of the time, one of those 
     things will be more important than any of the others, and contribute 
     disproportionately to the success of the design.

      - I also think it's useful to think of problems in terms of: "What's the 
        most important thing we're not getting right?"

      - That way, you don't need to focus on everything, just the most 
        important thing that's not being treated properly.

   - My current thinking is that H-bonds are that *thing* for DNA interfaces.
     
      - Almost all direct protein-DNA interactions are H-bonds.

      - H-bonds are notoriously difficult to design, due to rigid and narrow 
        geometries.

   - But, as Jacob argues, solvation and maintaining catalysis could both be 
     more important.

   - Solvation:
      
      - It may be that proteins just don't have enough freedom to position 
        their sidechain atoms to bind DNA without water.  In that case, trying 
        to design a dry interface may not be feasible.

      - People have put a lot of thought into water models, as Jacob has 
        described above.  I'd have to familiarize myself with that literature.

      - One idea I had is that you could maybe work out in advance the first 
        few solvation shells of DNA, and include the waters in those shells as 
        fairly static atoms in the simulation (maybe with partial occupancies).  
         
         - The more ordered a water is, the more important it is to simulate 
           explicitly.

         - I'm assuming the favored water positions for DNA are known.
     
   - Catalysis

      - I definitely agree that not perturbing recombinase function could be a 
        significant challenge.

         - This is also something that my library-vs-library approach will tell 
           me nothing about.

         - But we have a separate library-vs-library recombinase assay we can 
           use.  The in vitro assay should be better for getting directly at 
           the question of DNA binding, but we can switch to the recombinase 
           assay once we have some algorithms we're confident in.

      - I don't understand nearly enough about how recombinases work, so I'll 
        have to learn more about that before I start working on this is 
        earnest.

      - Machine learning might be useful here.

         - Surge has models that can learn what a "protein" is.
           
         - I wonder if we could train this model to learn what a "recombinase" 
           is, and then use it to identify mutations that aren't consistent 
           with being a recombinase.

         - In other words, I could do design normally and then use ML to filter 
           out designs that are predicted to interfere with catalysis or cause 
           other unintended large scale changes.
           
         - I imagine Surge would want to take this further and just use the 
           model to directly design recombinases, but that itself might be very 
           challenging.  There are certainly enough reasons why it's hard to 
           design recombinases, and ML doesn't really know about any of them.  
           Using ML as a filter might be an easier-but-still-useful 
           application.

- Learning from large-scale experiments

   - I agree that it's difficult to learn from failed designs, but I do think 
     it will be easier with more designs.  Each design might fail for its own 
     reasons, but with 10K designs, we should be able to see some patterns and 
     identify the most common failure modes.


2019/06/26:

Talk with Alex Garruss:

- The lab has had trouble getting these display techniques to work.

   - I should talk to Nancy; she apparently spent 2 years trying to get mRNA or 
     cDNA display to work.  Kettner should have her email address.

   - If I decide that I really want to get cDNA display to work, I should think 
     about trying to visit the lab that published those papers.  George can 
     afford it, and it'd be the fastest way to learn the technique.

- An alternative to B1H would be the following droplet-based assay:

   - Assay:

      - Construct DNA sequences with the following elements (in order):

         - DNA binding site
         - Promoter
         - Protein of interest (i.e. library member)

      - Express the designs in droplets using cell-free extract.

      - Hook up the protein of interest so that if it binds before the promoter, 
        it recruits the polymerase and results in the amplification of its own 
        sequence.

      - Sequence to see which designs are most represented.

   - The advantage of B1H is that it's an established assay.  This assay would 
     need to be developed.  

   - This assay would presumably be more quantitative than B1H, but not more 
     quantitative than the other in vitro assays I'm considering.

2019/08/02:

Talk with Kamesh:

- B1H is a good option:

   - I'm keeping B1H as a backup option, but Kamesh emphasized that it's not an 
     unattractive choice.  
     
   - One limit is transformation efficiency, but 1e7-1e8 could be enough 
     diversity to learn what I want to learn.
     
   - The assay is easy to get working.  I can order plasmids from Addgene, and 
     there are people in the lab (including Kamesh) who can help me gets 
     started.

- Consider benchmarking on existing datasets.

   - A lot of Zn-finger datasets are available.  It might be possible to create 
     a good benchmark.

   - That said, I really believe that the ability to test 1e4 designs will just 
     render benchmarks obsolete.

- Talk to people:

   - Martha Bulyk

      - Ex Church lab
      - Might have some ideas for applications besides recombinases
      - Might have some students that would be interested in helping.

   - Gary Storno
      
      - Friendly
      - Has done B1H

   - Scott Wolfe

   - Hamed Shateri

      - See Nature 2015/16 paper.

   - Tim Hughes

2019/11/18:

Jorge Marchand:

- To purify repA-DNA complex, try labeling the DNA with biotin.  The DNA is 
  likely to be available even if protein tags aren't, especially if the biotin 
  is on the 5' end (although this would require me to move Cy5 to the 3' end).

Even if protein tags are inaccessible, the Try labeling 

Ben Weinberg:

- B1H: Instead of transcribing GFP and doing FACS, transcribe a barcode and 
  read it using RNA-seq.

- And it's probably a good idea to just do B1H to get some data to work with.

2019/11/26:

Martha Bulyk:

2019/11/27:

- One potential application of protein/DNA interface design is to study 
  transcription factors.  Most transcription factors target many sequences, and 
  it would be nice to be able to deconvolute the roles of those sequences.  One 
  way to do that is to engineer TF variants that target a subset of the 
  wildtype TF's sites.  Similarly, this could be used to make synthetic TFs 
  with better defined specificities.

George Church:

- Nobody has really tried to fairly compare the off-target activity of CRISPR, 
  ZFNs, TALENs, and recombinases.  This would be valuable, especially for a 
  project trying to present itself as an alternative to CRISPR.  It isn't a 
  trivial comparison, either, because activity depends on the site being 
  targeted, the concentration of the protein, and the duration of editing.  The 
  location can't even be completely standardized, because all of these systems 
  recognize and act at sequences with slightly different offsets.

  But you could imagine measuring off-target activity for a number of 
  concentrations, times, and locations for each system, and comparing the 
  minimums.

- For my F32 proposal, I'll have to say why CRISPR won't just solve everything 
  (this echoes a comment made by Tanja).  Unfortunately, that means taking 
  potshots at CRISPR, and it's likely that some of the reviewers will have a 
  vested interest in CRISPR (since most people do, nowadays).  Prime editing is 
  probably my biggest competitor.  It's downsides are its efficiency, it 
  inability to make large insertions (it can insert loxP, but that requires two 
  steps to make a significant insertion), and its size (although the size could 
  be mitigated with better delivery options).  It's good to focus on the 
  efficiency of Cre, because that really stands out and doesn't rag on CRISPR 
  directly.

- George likes the chip-based approach to testing design algorithms.  I 
  mentioned that I was surprised by how slow this was catching on in the 
  protein design community, and he said he wasn't surprised because the 
  scientific community as a whole has been slow to adopt multiplexed DNA 
  synthesis since it was published in 2004.  Certainly the Church lab (and its 
  alumni) have been making good use of it, though.

- I might want to add Keith Joeng as co-collaborator, since he developed 
  GUIDE-seq.  But I'm now downplaying the role of GUIDE-seq relative to the 
  aims I sent George, and it's too late to do this kind of thing for the first 
  submission.  So if/when I have to resubmit, I'll ask Martha and Keith to 
  co-mentor.
  
