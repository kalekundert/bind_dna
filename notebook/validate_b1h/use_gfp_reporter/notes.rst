****************
Use GFP reporter
****************

[Noyes2008]_ (and every other B1H protocol I'm aware of) uses the HIS3/URA3 
reporters for auxotrophy-based survival screens.  I want to instead try using a 
fluorescent reporter.  This would have a few advantages:

- Less pressure to evade assay.
- Ability to directly observe non-binders (instead of inferring their presence 
  from what drops out of the library).

And disadvantages:

- Smaller libraries (maybe).
- Less coverage of hits.

Considerations
==============
- Gibson assembly is probably the way to go.

- Use p194 as template, and get rid of f1 ori while I'm at it.

- Get rid of f1 ori

- p27 has mWasabi, with Strep Tag.  p43 and p44 might have untagged mWasabi.  
  The original gBlocks I ordered with the fluorescent proteins don't seem to be 
  in my database.

- The different plasmid architecture I have require different Gibson primers.  
  Fortunately, I think it works out that I only need two primers:

  =========   ===========   =======   ======
  Construct   Zif268 Site   Barcode   Primer
  =========   ===========   =======   ======
  p193        in f1         yes       SR045
  p194        after f1      yes       SR045
  p186        before AmpR   yes       f1 5'
  p195        in f1         no        SR045
  p196        after f1      no        f1 5'
  p189        before AmpR   no        f1 5'
  =========   ===========   =======   ======

  - p193 and p194 have a terminator insert after the reporters.  There are two 
    annotated features between the reporters and the terminator: the f1 ori and 
    one of the rrnB terminators.  I want to remove both of these features 
    anyways: the f1 ori is irrelevant; the rrnB T1 terminator causes problems 
    for cloning (because it's a duplicate sequence) and is redundant since it's 
    followed by more terminators.
    
  - p186 and p189 have the Zif268 insert before AmpR, where it will have no 
    effect on cloning reporters.  For these constructs I want to keep the rrnB 
    T1 terminator, because otherwise the reporter won't have a terminator.  But 
    I can still take the opportunity to delete the f1 ori.  The primer will 
    have to overlap with the f1 ori a little bit, though, because the 
    terminator doesn't make for a good primer.

    Note also that I think the RepA gene lacks a terminator in these 
    constructs, which could be a problem.
    
  - p195 doesn't have the barcode/terminator insert that p193 and p194 do, but 
    it does have the same SR045 cloning scar.  This will put the reporter 
    back-to-back (although facing the opposite direction) as the Zif268 gene, 
    with no terminators in either direction.  I'd be uneasy using that 
    construct for real, but that's what p195 is.

  - p196 still has the rrnB T1 terminator between the reporter and the Zif268 
    gene, so I can use it in the same way I'm planning to for p186 and p189.  
    It also has the SR045 cloning scar, but using it would remove the rrnB T1 
    terminator, which I think would be a bad idea.

Internal control
----------------
One big advantage of not using a survival reporter is that I can include a 
internal control.  With GFP as the primary reporter, RFP would be the natural 
internal control.

- The RFP would point in to opposite direction as GFP, such that GFP would 
  remain the only gene on the plasmid facing its direction.

- I'm going to use mCherry.  [Piatkevich2011]_ recommends mRuby, mCherry, and 
  mStrawberry as the most prominent RFP variants.  mRuby is the foremost 
  recommendation, but the other two have ≈10x faster maturation times.

  - The strain I used for my sgRNA project had mRFP1 and sfGFP.

- Instead of doing my own reverse translation, I searched Addgene for bacterial 
  expression plasmids with mCherry.  The top hit was #29769.  I plugged it's 
  mCherry sequence into a codon usage tool 
  (https://www.bioinformatics.org/sms2/codon_usage.html).  It doesn't look 
  like a very sophisticated reverse translation: it pretty much just uses one 
  codon for each amino acid.  But presumably the sequence works and expresses 
  well, so I won't worry about it.

  - This mCherry sequence does have a BsrGI site, though, so I'll have to get 
    rid of that.

- The Salis Lab has a reverse translation tool (called CDS calculator) that 
  seems good.  It accounts for codon preference, restriction enzymes, 
  repetitive sequences, internal RBS, internal promoters, internal terminators, 
  and other things like that.  Regarding restriction enzymes, I chose to 
  exclude the enzymes I uses (or would like to use) to change the target site: 
  EcoRI, BsrGI, and NotI. Here's the sequence it predicted for mCherry:


Promoters/RBSs
--------------
GFP will have the same promoter/5'-UTR as the HIS/URA reporter, of course.

For RFP:

- I could use the same promoter/5'-UTR as GFP:
  
  - That would give equivalent numbers of GFP/RFP molecules in the case of no 
    binding.  But I can't think of why that would be useful.  The brightness 
    levels are different.

  - I don't necessarily want RFP expression to be that low; it could be hard to 
    measure reliably.

- I could use the constitutive promoter from [SegallShapiro2018]_: J23105

  - These promoters were used with pSC101 and with a fluorescent protein, and 
    were shown to be pretty proportional to copy number.

  - Sequence: GGCGCGCCTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCAAGGT 

  - J23105 is a medium-strength promoter, according to IGEM.

  - This promoter includes an RBS too: B0032

    - Sequence: TCACACAGGAAAGTACTAG

    - This also seems to be a medium-strength RBS.

    - Note that the sequence on the IGEM wiki for this part is missing the 6 3' 
      base pairs.  I trust the sequence from the paper though, because I can 
      see that it was actually used in these constructs.

  - It wasn't totally clear to me from the SI whether there were any spacer 
    sequences between the indicated parts.  But I confirmed that there are no 
    spacers (i.e. the parts go together back-to-back) by looking at pTHSSe_31, 
    which has a complete sequence on addgene.

- I could use promoters from the Salis Lab non-redundant parts list.

  - These include sgRNA target sites, apparently so that Cas9 can regulate 
    expression somehow.  I'm not familiar with that technique, and I certainly 
    won't use it for this project..

Should I use the lac operator?  The other genes all have it, but I can't really 
understand why.

Terminators
-----------
- The HIS/URA reporter uses a partial rrnB terminator.  This causes issues with 
  cloning, because the complete rrnB terminator is used upstream of the 
  reporter.  

- I think I should just use the [Chen2013]_ terminators.  I've already been 
  using them, and I just realized that they're the terminators used by the 
  Salis lab (even though they come from the Voigt lab), so I think they're 
  regarded as good options.

  - I've already used:

    - L3S2P21
    - ECK120033737
    - ECK120029600

  - These are the 3 strongest terminators in the most recombination-resistant 
    set published by [Chen2013]_.  See Supplementary Tables 2 and 3.  If I need 
    more terminators, I don't see any reason to not just keep taking from that 
    set.

RiboJ
-----
Reading these synthetic biology papers, I see that it's common to include RiboJ 
between the promoter and the RBS.  I didn't know what RiboJ was, so I had to 
research it.  [Clifton2018]_ seems like a good introduction.  It's meant to 
give transcripts a consistent 5' end regardless of which promoter is used, but 
it seems to increase gene expression (possibly by increasing the stability of 
the mRNA or by helping keep the RBS exposed).

If I were to use this, I'd only use it for one gene, because I don't want any 
duplicate sequences.  It also seems to make the most sense when long, natural 
5'-UTRs with possibly cryptic functions are being used.  It seems less useful 
in a very synthetic situation like this.
