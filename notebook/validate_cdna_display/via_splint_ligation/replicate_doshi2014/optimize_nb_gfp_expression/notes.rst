**************************
Optimize Nb-GFP expression
**************************

My goals for this experiment are:

1. Show that Nb-GFP is expressed well in PURExpress.
2. Determine the optimal DNA and mRNA template concentrations.

.. protocol:: 20220103_make.txt 20220103_optimize_nb_gfp_expression.txt

.. figure:: 20220103_titrate_f145.svg

Observations:

- The gel is smiling a lot.  That didn't happen in :expt:`117`; not sure what 
  the problem is.

- The RNA lanes have bands at ≈8 kDa.

  - I plugged the mRNA sequence into the Salis Lab RBS calculator, to see if 
    it's possible that this is an internal-RBS product.

  - The main RBS has a predicted translation rate of 1e5.  The next highest 
    site is Met63, with a predicted rate of 1e2.  That's not very high, but it 
    does have about the observed MW (8 kDa).

  - I can't really tell if these bands are present in the negative control.

- I can't really say whether or not Nb-GFP was expressed.  If I had to guess 
  I'd say that the 14-17 kDa bands are somehow indicative of Nb-GFP expression, 
  but I'm not sure.

Next steps:

- I should try western blotting.  My Nb-GFP constructs are FLAG tagged, so I 
  can use the antibodies I already have.


