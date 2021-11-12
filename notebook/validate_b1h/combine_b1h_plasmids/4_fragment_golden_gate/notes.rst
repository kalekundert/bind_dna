**********************
4-fragment Golden Gate
**********************

2021/09/10:

I want to try cloning it via a series of conservative steps.  This should 
reduce the chance of ending up with a defective plasmid, and help determine 
what the problem is if I do.

Here are the steps I envision:

- Replace KanR with AmpR.  This isn't strictly necessary, but I strongly prefer 
  working with Carb.

- Install Zif268-rpoZ

- Remove f1

- Remove unannotated regions.

Considerations
==============

Where to insert Zif268-rpoZ?
----------------------------
I want a site that is as unlikely as possible to cause problems.  Candidates:

- Between URA3 and f1 ori

  - Don't really have a primer that can do this.

- Inside f1 ori

  - I can't really know if unannotated sequences are important or not, because 
    I don't know if they do anything.  But I know what the f1 ori is supposed 
    to do (enable viral packaging for mammalian expression) and I know it's not 
    important to me.  That's not to say that it couldn't have some cryptic 
    secondary function, but it's a good candidate "landing pad".

- Between f1 ori and Rep101

  - Unannotated
  - Leaves the HIS3/URA3 mRNA unchanged.

- Next to the target

  - May not be a good idea for reasons discussed in :expt:`-1`, but worth 
    trying for the sake of completeness.

Results
=======

.. protocol:: 20210923_pcr_verification.txt

.. figure:: 20210923_check_p184_p185_p186.svg

- The negative control for the o188 reaction wasn't clean, so I should be 
  skeptical of those results.

- None of the picks seem to be correct.

- I'm not sure why most of the reactions have a smear around 0.2 kb.

- I think this is a pretty good test.  I'll try the cloning again with 
  gel-purified fragments and more competent cells.

.. protocol:: 20210928_check_rpoZ_zif268_insert.pdf 20210927_make.txt 20210928_check_rpoZ_zif268_insert.txt

.. figure:: 20210928_check_rpoz_zif268_insert.svg

- Got lots of p184 colonies; very few p185/p186 colonies.

- Still no clear hits.

2021/10/13:

After checking the hits I'd found previously, and finding that none had the 
Zif-rpoZ insert, I decided to start screening colonies using o266:

.. protocol:: 20211013_check_rpoZ_zif268_insert.txt

.. figure:: 20211013_check_p184.svg

- Still no hits, although I'm tempted to sequence the plasmid with the 
  insert...

2021/10/21:

.. protocol:: 20211021_check_rpoZ_zif268_insert.pdf 20211021_check_rpoZ_zif268_insert.txt

.. figure:: 20211021_check_p184.svg

- The bands are hard to see, but there are a few hits with ≈1kb amplicons.  
  None with the expected insert (2.5 kb amplicon), though.

- I'm not sure what else I can do to improve the assembly.  I already used 
  gel-purified fragments and lots of competent cells.  Ideas:

  - No terminators: just a single insert.

  - Clone GFP instead of zif: I can see it, and I can use it as a drop-out in 
    subsequent steps?  Seems like a lot of work, though.

