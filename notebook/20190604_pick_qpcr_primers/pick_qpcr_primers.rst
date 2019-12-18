************************
DNase: Pick qPCR primers
************************

.. update:: 2019/12/18

   I redesigned my plasmids to explicitly include [Subramanian2018]_ PCR 
   primers, so the results of this experiment are no longer relevant.

I am planning to use qPCR to measure whether or not Zif268-binding can protect 
a DNA barcode from DNase treatment.  Ultimately I will want to make this 
measurement using NGS, but qPCR makes more sense for the initial, 
low-throughput, control experiments.

The first step is to find a set of primers and primer concentrations that work 
well for qPCR.  The plasmids I'm working with initially (11 and 15) include an 
explicit reverse primer (skpp-202-R), but no forward primer.  This was an 
oversight on my part, and maybe something I'll correct later.  For now, I just 
designed forward primers to anneal with various parts of the target and barcode 
sequence.

2019/06/04:

I designed primers manually such that each:

- is 20 bp long
- is between 40-60% GC
- has a 3' G or C

The primers anneal to the "Zif268 target" and "DNA barcode" features of the 
plasmid.  The primers are named for the length of the resulting amplicon (e.g.  
"QPCR_51_REV" is a smaller amplicon than "QPCR_57_REV").

.. datatable:: primers.xlsx

2019/06/06:

I did standard PCR to check for clean amplification.  I decided against doing a 
temperature gradient, because I'm not sure that Q5 and the Biorad polymerase 
prefer the same temperatures.  I can do a temperature gradient later on the 
actual qPCR machine.

.. protocol:: 20190604_pcr.txt

   - I used pDBP011 as the template.

.. figure:: 20190606_validate_qpcr_primers.svg

   neg: The negative control, which included template but no primers.  51, 57, 
   60, 67: The primer pairs.  Note that the numbers represent the length of the 
   expected product for each primer pair.  The lanes left of the ladder are 
   from another experiment.

- I wish I'd included a "primer only" lane, to be sure that I'm really seeing a 
  short amplicon and not just unreacted primers.  But I do think size of the 
  bands (as measured by the ladder) is more consistent with a short amplicon 
  than with unreacted primers.

- All of the reactions seemed to be clean and efficient.  But I only want to 
  move forward with one, so I used ImageJ to quantify which lane had the most 
  product.  On this basis, I chose the "60" primer pair.  This could be unfair, 
  because the different intensities could be due to the reactions having 
  different amounts of primer (the primer volume was small and I wasn't careful 
  about pipetting accurately).  But it doesn't really matter if this isn't the 
  optimal choice, as long as it works well.

.. protocol::

   - Subtract background (50 px rolling ball)

   - Draw lanes

   - Count pixels

.. datatable:: 20190606_validate_qpcr_primers.xlsx

   Intensity of the amplified product band for each primer pair, as measured by 
   ImageJ in units of pixels.

