*********************
Replicate [Doshi2014]
*********************

Erkin has successfully coupled protein and mRNA using the protocol described by 
[Doshi2014]_.  Given this, I'm very hopeful that I can get the same protocol to 
work.  This protocol is very similar to the others I've tried, with a few 
differences:

- Splint ligation.
- λ exonuclease is added to digest the splint and any unreacted linker.
- Magnetic oligo-dT(25) beads are used to purify the mRNA after both ligation 
  and translation.

I'm going to start by trying to reproduce the protocol from the paper exactly, 
i.e. using all the same reagents.  Erkin has these reagents, so I don't need to 
order them myself.  Once I can get that to work, I'll have to figure out how to 
attach my RT primer (and maybe some fluorescent moieties).

Notes
=====

Splint modifications
--------------------
Here is the splint used by [Doshi2014]_::

  /5Phos/TTTTTTTTTTTTTTTAGAACCACCACCAGAAC/3ddC/

I was initially confused as to the purpose of the 5' and 3' modifications.  But 
I figured it out:

- 5' phosphate: This makes the linker a substrate for λ exonuclease.  According 
  the NEB: "Preferred substrate [for λ exonuclease] is 5'-phosphorylated 
  double-stranded DNA although non-phosphorylated substrates are degraded at a 
  greatly reduced rate. [...]  λ exonuclease is ideal for conversion of linear 
  double-stranded DNA to single-stranded DNA via preferred activity on 
  5'-phosphorylated ends."

  In principle, this modification also makes it possible for the splint to be 
  ligated to the mRNA, but my guess is that this just doesn't happen because 
  there's nothing to hold it in place.

- 3' ddC: This is to prevent linker/splint from being ligated to the splint.  
  Since neither the linker (3' puromycin) nor the splint (3' ddC) has a 3'-OH, 
  the mRNA is the only possible 5' ligation partner.  It's particularly 
  important that the splint and linker don't react, since both are in 5x excess 
  (IIRC).

Splint mismatch
---------------
In contrast to [Keefe2001]_ and [Barendt2013]_, the splint does not perfectly 
align to the ends of the mRNA and the linker:

- [Doshi2014]_::

    mRNA                                         puro
    UGUCGUCCGGUGGUUCUGGUGGUGGUUCUGGUGGUGGUUCUGGU GCAAAAAAAAAAAAAAAAAAA/iSp18//iSp18/ACC/3Puro/
                       /ddC/GTTCTGGTGGTGGTTCT------AAAAAAAAAAAAAAA/Phos/
                            splint (reverse complement)

- [Barendt2013]_::

    mRNA (pRDV_BbsI_r)   puro
    GCTTGGGTCTTCCGGTGGCC AAAAAAAAAAAAAAAAAAAAAAAAAAACC/3Puro/
             TTCCGGTGGCC-AAAAAAAAAA
             splint (reverse complement)

- [Keefe2001]_:
    
  This is more of a review, but calls for a splint with 10nt of complementarity 
  to both the mRNA and the puromycin linker (which is almost exactly what 
  [Barendt2013] uses).  

I don't know if there is a reason for this.

