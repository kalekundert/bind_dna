*********************
Validate cDNA display
*********************
I want to see whether or not I can successfully attach a translated protein to 
the cDNA encoding it.  This is the first step in doing cDNA-display-based 
DNA-binding assays.  My goal is to run a gel where both the protein and the 
cDNA can be seen in different channels, and to show that the two bands are 
superimposed.  As a negative control, I can express the protein without the 
linker.

Considerations
==============

Visualization
-------------
To see whether or not CIS-display was working, I expressed GFP using DNA 
labeled with Cy5, so I could simultaneously visualize both molecules.  See 
:expt:`36` for details.  To see whether or not cDNA display is working, I'll 
have to do things a little differently, because the linker I ordered is labeled 
with FITC.  This means that I have to visualize the protein in the Cy5 channel.  
I thought of a few ways to do this.  Of these ideas, I think the fluorescent 
protein stain seems the most promising:

- Use a fluorescent protein stain.
  
   - Because cDNA display creates a covalent attachment, I should be able to 
     use SDS-PAGE to see whether or not the attachment was formed.

   - Because SDS-PAGE is so much higher resolution than native PAGE, I will 
     probably be able to see a distinct band for my product even if I just 
     stain all proteins.  This would remove the need to do any purification 
     steps (although such steps would still give nicer data).

   - Because I wouldn't have to use a fluorescent protein, I could just use my 
     Zn-finger, which would make my result more relevant to the actual binding 
     assay I want to perform.

   - While I could use Coomassie, a fluorescent stain will be more sensitive, 
     and will allow me to create nice multichannel images using the gel 
     imager.
     
   - A brochure for the :download:`Sapphire imager <sapphire_brochure.pdf>` 
     claims that SYPRO-Ruby can be visualized using the Cy5 laser.  I'm a 
     little skeptical because SYPRO-Ruby doesn't appear to be excited by 658 nm 
     light, but I trust the brochure enough to give it a try.
     
   - SYPRO-Ruby might overlap with the GFP channel.  It will be almost 
     maximally excited by the 488 nm laser, but it shouldn't emit in the 
     507-529 nm band.  If there is cross-talk, I might also try SYPRO-Ruby, or 
     maybe the Typhoon will have better laser/filter combinations.

- Use a near-IR fluorescent protein:

   - I don't want to use an RFP, because most RFPs have some overlap with the 
     GFP channel.  Even if the overlap is small, it makes the results very 
     difficult to interpret.

   - Near-IR proteins all seem to require the biliverdin cofactor.  This could 
     cause issues with background fluorescence or getting the protein to fold, 
     although most likely it would be fine.

   - I would need to continue using native gels, so that the fluorescent 
     protein would remain folded.

- Western blot

   - Kinda the best of both worlds, in that I can run an SDS-PAGE gel and then 
     visualize just my protein of interest.

   - I'd probably just use an epitope tag to find my protein, which I could 
     possibly also use for purification (see below).

   - I haven't done a lot of Westerns before, so I'd have to learn how to do 
     that.  Could be a bit of a learning curve, but it'd also be a good 
     technique to know.

- Order a new linker with Cy5 instead of FITC.

   - This would be expensive.

   - But it would make it possible to visualize the mRNA/linker conjugation 
     using SYBR green II.

   - But then I couldn't use SYPRO-Ruby, which could be really nice.  I could 
     use SYPRO Orange, though, it shouldn't overlap with Cy5 at all.

     .. note:: 
     
         Previously I was planning to use LUCY-506, which is more of a green 
         stain.  But it is discontinued, and doesn't seem to be available 
         anywhere.

   - If I were to do this, though, I should also think about ordering the 5' 
     end phosphorlyated.  That would let me skip the phosphorylation step, and 
     would only improve efficiency.

   - I could do a test run by ordering an unbranched Y-tag oligo from IDT and 
     seeing if it ligates well.

Purification tags
-----------------
Since I'm planning to do SDS-PAGE, it shouldn't be absolutely necessary to 
purify my proteins.  But I am worried about getting partially transcribed 
products due to puromycin being present in the reaction.  To account for this, 
I'd like to have a C-terminal tag I can use to purify fully-translated 
products.  

Note that this C-terminal tag will really be an internal tag, because the 
puromycin linker will be attached to the C-terminus.  This limits the selection 
of tags I can use.  I consulted the list of affinity tags compiled by 
[Kimple2013]_ to find the best candidates.

- His-tag: [Naimudden2016]_ used a C-terminal His-tag (flanked by GGGS motifs 
  on either side) and seemed to be able to successfully purify protein, even 
  though [Kimple2013]_ claims that His-tag are not internal tags (although I 
  don't see why they would need to be at either terminus to coordinate Ni).  Of 
  course, all of the proteins in the PURExpress reaction are also His-tagged, 
  but I could always use a 100 kDa spin-filter get rid of most of them 
  (although not any that might run near my protein in a gel).

- Strep-tag: I don't really have any reason to believe that this would work 
  internally, except that the Twin Strep-tag has one Strep-tag that's basically 
  internal, and it still seems to contribute something to binding.  Plus I 
  already have all the beads and buffers, so I might as well try it.

- HA-tag, Myc-tag: These are two common epitope tags that are considered by 
  [Kimple2013]_ to be internal tags.  Antibody purifications are considered to 
  give low yields, but I think they're worth trying if the above tags don't 
  work.

