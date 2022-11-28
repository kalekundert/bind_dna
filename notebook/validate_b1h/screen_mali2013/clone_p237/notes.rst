**********
Clone p237
**********

2022/09/06
==========
.. protocol:: 20220906_make_p237.pdf 20220906_make_p237.txt

  - f186: don't have gel-purification protocol yet, and haven't used my 
    electroelution devices yet.

  - f187: digest protocol doesn't work well with small volumes, need to check 
    that by hand.

  - l3: Current protocol is to do 10 cycles, then more until a band is visible 
    (i.e. ≈1 ng/µL).  Fitzy recommends qPCR, and that does sound more 
    quantitative to me.  But I'll start with this and see what happens.

.. datatable:: 20220909_electrotransform_p237.csv

- I only got 1e6 transformants.

  - I only expect this library to have 1e4 members, so 1e6 transformants will 
    be enough.

  - 1e6 is not very good, though.  There's definitely room for improvement, and 
    that improvements might be necessary for future steps.

  - Ideas to get more transformants:

    - Gel-purify the insert.
    - Use more backbone/insert.
    - Use more cells.
    - Try home-made cells.

2022/10/18
==========
After sequencing f202 and f203, I'm not confident that any of the libraries 
I've cloned are correct, because I didn't use phosphatase or perform −insert 
control transformations.  I'm going to sequence everything I have, to make sure 
I'm not operating on any false assumptions, and in the meantime I'm going to 
get ready to remake everything.

2022/10/19
==========
The sequencing reactions from yesterday all failed.  While I don't necessarily 
think that means the constructs are wrong, I'm going to remake them anyway.

In the future, one reason to prefer Nanopore sequencing is that it will 
probably be more informative for totally wrong constructs.

.. protocol:: 20221019_make_f186.txt

.. figure:: 20221019_gel_purify_f186_f190.svg

- I purified the bottom part of the bottom band.  My recovery wasn't very good, 
  but I was trying to be selective.

- I'm not sure what the smeary high-MW bands are.

  - The whole plasmid is 2.4 kb, so I can't explain the 6 kb band.

  - It may be that digestion did not go to completion.  Nicked DNA runs slower 
    than linear DNA, so that could explain the 6 kb band.  If this were the 
    case, though, I'd expect to see a supercoiled band running faster than the 
    linear one.  So I'm not sure that this is what's happening.

  - The only other alternative I can think of is something is wrong with p236.

  - See :expt:`207` for follow-up on this.

2022/10/24
==========
.. protocol:: 20221024_make_p237.pdf 20221024_make_l3.txt 20221024_make_f187.txt 20221020_make_p237.txt

.. figure:: 20221025_electrotransform_p237.svg

Observations:

- I think that most of the transformants just have undigested p236 (not p237).

  - I got the same number of colonies with/without insert (f187), which means 
    that either (i) I have undigested backbone or (ii) the backbone is 
    re-ligating.

  - I don't think that re-ligation is the problem, because I added phosphatase 
    this time.

  - I see from the 10/21 gel that uncut p236 runs almost identically to f187.  
    I'm suspicious of this result (for reasons described above), but if I take 
    it at face-value, it explains why the gel purification may not have done 
    anything.

Next steps:

- Make f186 by PCR.
  
  - This would allow me to use a very small amount of p236 to start with, and 
    the digest what I do use with DpnI.  Presumably, this would pretty much 
    eliminate the possibility of transforming undigested p236.

  - I'll keep the gel purification step, to guard against unexpected 
    amplification products.

  - I'm tempted to use fluorescent primers to monitor the restriction 
    digestion, but that's probably overkill for now.

2022/10/28
==========
.. protocol:: 20221026_make_f186.pdf 20221027_make_p237.pdf 20221026_make.txt 20221027_make.txt

  - Midiprep:

    - Resuspended cells in 2x5 mL water.  Recovered ≈5 mL total.

    - Measure A600(10mm)/100: 0.41

    - Assume A600(10mm) of 0.1 implies 1×10⁸ cells/mL.

    - Calculated that 50 mL overnight culture (assumed density: 3-4×10⁹ 
      cells/mL) would correspond to 4.87 mL of resuspended cells.  Used all 5 
      mL.

.. figure:: 20221028_electrotransform_p237.svg

.. datatable:: 20221028_electrotransform_p237.csv

- I should probably keep using 0.5x TA buffer for gel purifications:

  - I used LAB buffer when gel-purifying f186 this time.
  - I used 0.5x TA when gel-purifying f186 on 10/19.
  - The TA gel had sharper and more well-resolved ladder bands.
  - Perhaps I could find parameters that would make LAB work better.  For 
    example, I used a much higher voltage/shorter time for LAB than I did for 
    TA.  But the TA results aren't bad, and it's probably not worth the hassle 
    to optimize LAB.

- Using PCR to amplify the backbone gave a very low level of background

  - This adds more support the idea that problem was getting rid of undigested 
    plasmid, not re-ligation.  I already was pretty confident in this idea, see 
    the 10/24 data.

  - Of course, using PCR will not be an option for the next round of cloning 
    steps (with BsmBI).  Neither will be using a large insert in order to 
    facilitate gel purification.  I'll have to get the digestion to go to 
    completion.

2022/11/02
==========
I sent the plasmid I purified on 10/28 for full-plasmid sequencing, but the 
reaction failed.  Specifically, there were only 1-2 reads.  I don't think the 
problem was that I didn't send enough plasmid, because I measured concentration 
via Qubit.  However, there was ≈100 ng/µL of RNA in the sample.

Today I treated the plasmid with RNase, which appeared to successfully remove 
all the RNA.  (The Nanodrop and the Qubit now give nearly identical readings.)  
I sent the plasmid for sequencing again, this time sending two concentrations 
(30 ng/µL and 300 ng/µL) to cover my bases.

I also did a handful of test digests to see if the plasmid seems right:

.. protocol:: 20221102_test_digest.pdf 20221102_test_digest.txt
.. figure:: 20221102_test_digest_p237_p239.svg

- Maybe I should've done double-digests, to get different banding patterns for 
  the different enzymes.  Of course, the digestion don't all go to completion, 
  and an incomplete double-digest would be harder to interpret.

- p237 has the expected banding pattern.  Mostly supercoiled for the negative 
  control and EcoRV, 2.4 kb for everything else.

- Based on this, I'm going to assume that the sequence is right and move ahead 
  with making p241.

2022/11/03
==========
I got the Sanger sequencing results, and I realized that I ordered an old 
version of l3, different from what I have in my SnapGene files.  The older 
version has an AGAT BsmBI junction, along with half of sr091 (to be used when 
amplifying the barcode).  The current version instead has an ATGG junction and 
no extra base pairs.  These base pairs make up the end of sr116, and as a 
result this approach results in a shorter overall library.

The easiest way to paper over this would be to make "f196" by PCR, with primers 
that install a different BsmBI overhang.  I'd end up using different primers to 
amplify the barcode in the end, but that's not a big deal.

Alternatively (or additionally), I could reorder the library.  That would take 
a few days to arrive, but it would let me add more library members (since I 
regret choosing just 4 sites to test).

Discussion
==========
So far, I've had two major problems with this experiment:

- Plasmid preps give very low yield, see :expt:`194`

  - Using Zymo kits seems to help.

- Restriction enzymes have very low activity, see :expt:`207`.

  - Doing long incubations with a large excess of enzyme seems to help.
