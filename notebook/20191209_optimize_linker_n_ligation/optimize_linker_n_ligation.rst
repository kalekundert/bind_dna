**************************
Optimize linker-N ligation
**************************

[Naimudden2016]_ reports a 90--95% efficiency ligating linker-N to mRNA.  I 
want to achieve similar efficiencies.

Results
=======

Ligate linker-N --- 2019/12/9
-----------------------------
.. protocol:: 20191209_transcribe_rna.txt

   See binder: 2019/12/9 and 2019/12/13

.. figure:: 20191213_ligate_linker_n.svg

- The transcribed RNA is not very homogeneous.  See 
  :expt:`20191216_optimize_mrna_transcription` for more discussion.

- The ligation was 66% efficient, less than the 90--95% efficiency reported by 
  [Naimudden2016]_.  But I have a number of things I can try (discussed in the 
  :ref:`validate_cdna_display_ligation` section) to improve this.

  Note that this efficiency is probably a slight overestimate.  I calculated 
  efficiency using the same equation as [Naimudden2016]_, but this equation 
  doesn't account for the fact that the conjugate has 28 bp of double-stranded 
  DNA/RNA hybrid.  `According to Biotium 
  <https://biotium.com/faqs/gelred-gelgreen-ssdna-rna/>`, "titration assays 
  using a fluorescence microplate reader showed that the fluorescence signal of 
  GelRed® bound to ssDNA and RNA is about half that of GelRed® bound to dsDNA."  
  Assuming that double-stranded DNA/RNA is as bright as dsDNA, this would give 
  a corrected efficiency of 64%.

  There are also reasons why this efficiency could be just plain inaccurate.  
  One is that the smeary RNA made subtracting the background rather subjective.  
  Hopefully I can improve this by getting cleaner RNA.  Another is that there 
  could be some FITC signal in the red channel.  To check for this, I need to 
  measure both the red and green channels before adding GelRed, which I didn't 
  do this time.  Note that the efficiency looks much lower in the 300 nm GelRed 
  image.  This image shouldn't have any signal from FITC (another thing I 
  should test), but it does have a smear that could be making the lower band 
  seem brighter.

- Next time I do this experiment, I should setup control reactions without 
  linker and mRNA.  This way, all three lanes would have the same amount of 
  material, which would make the gel easier to interpret.

- Linker-N runs about with the dye front.  So don't run the dye front off the 
  gel next time.  That said, I'm mostly interested in the difference between 
  the two mRNA bands, and running the gel longer might help resolve them 
  better.

- Note sure what that high-MW linker-N band is.  (It's more easily seen in the 
  "intensity level 3" image that I didn't include here.)  But it might be a 
  consequence of the lane being severely overloaded.

- I think the green scratch is caused by the EZdoc UV tray.  The laser scanner 
  images without the scratch (not shown here) were taken before I'd added 
  GelRed or imaged with the EZdoc, and the image with the scratch was taken 
  after.  I thought the scratch could also be due to something on the bottom of 
  the tip-box scratching the gel during shaking, but the scratch (vertically 
  all the way from top to bottom, rather than circular) is not really 
  consistent with that.  

Ligate Linker-N
---------------
- Image the gel using the 488 nm and 520 nm lasers (Sapphire) and the 300 nm 
  illuminator (EZdoc) before adding any GelRed.  This will allow me to be more 
  confident about overlapping signal between fluorescent channels.
