*************************
Validate ligation via Cy5
*************************

[Doshi2014] doesn't include any fluorescent labels in the puromycin linker.  
This makes it hard to tell it the ligation reaction is really working, because 
the difference between the ligated/unligated bands is very small.  It's also 
hard to tell the difference between annealed and ligated oligos.

All of these problems would be solved if the puromycin linker was Cy5-labeled, 
because then it's presence could be detected regardless of how much it causes a 
band to shift.  Of course, it's possible that Cy5 itself could interfere with 
the ligation (or more likely, the downstream coupling step).

Considerations
==============
- The [Doshi2014]_ linker::

    /5Phos/GCAAAAAAAAAAAAAAAAAAA/iSp18//iSp18/ACC/3Puro/

- [Reyes2021]_ includes FITC in the puromycin linker::

    /5Phos/CCCTTCACCTGATCCGCTGAAAAAAAAAAAAAAAAAA/iSp18//iSp18//iFluorT//iSp18/CC/3Puro/

- The [Nagumo2015]_ linker::

    p(dCp)2-T(fluorescein)p-PEGp-(dCp)2-puromycin

  This protocol doesn't appear to use splint- or Y-ligation.  Just mix the 
  (6-mer) linker and the mRNA with RNA ligase.  I'm a bit skeptical of that.

- The [Naimudden2016]_ linker::

    CCCCCCCGCCGCCCCCCG(5-Me-dC-(5'C)(5'C)(5'T)(5'T)(5'G))A18(Spec18)(Spec18)(Spec18)(F-dT)(Spec18)CC(Puro)

There seems to be precedent for putting fluorescein-dT one spacer-18 before the 
3' CCA motif.  Cy5 has a fairly different chemistry from FITC (it's not 
attached to a nucleotide; it participates directly in the phosphoramidite 
reactions), but this seems like as reasonable a starting place as any::

  o297: /5Phos/GCAAAAAAAAAAAAAAAAAAA/iSp18//iCy5//iSp18/ACC/3Puro/




