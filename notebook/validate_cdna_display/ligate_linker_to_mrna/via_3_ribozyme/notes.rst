***************
Via 3' ribozyme
***************

:expt:`1` shows that only about half of the mRNA appears to be capable of being 
ligated to the linker.  This maybe due to heterogeneity at the 3' end of the 
mRNA.  There are several ways to eliminate that heterogeneity [Avis2012]_, but 
the simplest and most popular is to include a cis-acting ribozyme after the 
sequence of interest [Schurer2002]_ [Walker2003]_ [Akoopie2018]_.  The HDV 
ribozyme is efficient and puts no sequence requirements on the transcript.

Considerations
==============

Which ribozyme?
---------------
There are a couple different variants of the HDV ribozyme that have been 
described in the literature:

.. datatable:: hdv_ribozymes.xlsx

- [Schurer2002]_:

  - 67 nt
  - Probably the most widely used.
  - Least efficient if the nucleotide 5' of the cleavage site is an A, which it 
    is in my case.  But the effect isn't large (71% vs 79% cleavage).

- [Akoopie2018]_:

  - 56 nt
  - Reported to be better for some sequences and worse for others, but slightly 
    worse that [Schurer2002]_ on average.  Might be worth trying both.

I think I'll start with [Schurer2002]_.

Terminator
----------
Currently, I digest the plasmid encoding the mRNA to prevent RNAP from reading 
beyond the end of the Y-tag.  With a 3' ribozyme, this step may not be 
necessary.  Instead, I can add a T7 terminator after the ribozyme.  Termination 
won't be 100% effective, but it doesn't matter because every transcript will 
get processed to the right size anyways.  

While this would remove the need for the restriction digest, it might also 
introduce the need for a phenol-chloroform extraction (to remove RNase from 
miniprepped plasmid DNA).

.. update:: 2020/10/28

  In the interest of being able to clearly see if the ribozyme worked as 
  intended, I digested the plasmid with HaeII.  This should lead to the 
  creation of two fragments with defined length: the gene (829 bp) and the tail 
  (347, 336 bp for HDV67, HDV56).

Note that the whole terminator is transcribed, so I don't need to worry about 
including a spacer between the ribozyme and the terminator.  Briefly, 
rho-independent termination works like this [Peters2011]_:

- The terminator consists of a GC-rich hairpin followed by a U-rich tract.
- RNAP pauses on the U-rich tract.  I'm not totally sure why.
- This allows time for the GC-rich hairpin to form.
- The hairpin dislodges the RNAP.

Cloning
-------
Cloning in the Y-tag region is tricky, because there's a lot of GC- or AT-rich 
sequences.  I think the best bet is to order the HDV/T7-terminator cassette as 
a gBlock, then install it using Gibson assembly.  I decided to use XmnI to 
linearize the backbone.  This prevents the inserts from being specific to the 
affinity tag (e.g. His6, Strep, etc.) and avoids any difficult PCR steps.  It 
does keep the AAAA overhang created by XmnI, rather than the AAA overhang used 
by [Naimudden2016]_.  I don't think this overhang is a problem, though.  And in 
any case, it avoids changing two variables at once.
