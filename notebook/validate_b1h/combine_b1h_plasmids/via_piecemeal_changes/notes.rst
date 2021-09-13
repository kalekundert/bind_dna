*********************
Via piecemeal changes
*********************

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
