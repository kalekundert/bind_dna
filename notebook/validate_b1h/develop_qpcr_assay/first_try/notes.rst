*********
First try
*********

2022/05/05:

.. protocol:: 20220505_qpcr_assay.pdf 20220505_qpcr_assay.txt

.. figure:: 20220507_try_b1h_qpcr_assay_heatmap.svg

  Heatmap showing raw Cq values from the experiment.  See the layout file in 
  the experiment directory for a description of the contents of each well.

.. figure:: 20220507_try_b1h_qpcr_assay.svg

  Bar plot showing the relative expression of the 2TGG and 2AAA targets for the 
  4 plasmid architectures I built and tested.

Observations:

- There is a clear difference between target and non-target, but the difference 
  is only 4x.  I was expecting higher, based on the auxotrophy results.

- I didn't include the -RT controls, so I don't know if any of the signal is 
  due to the plasmid itself being amplified.

- I'm not sure if the differences between the plasmid architectures are 
  significant, but the best performing one (p224) is the same architecture that 
  performed best with the HIS/URA selection (p189).  This is also the construct 
  that presumably has the highest level of AmpR expression.  This assay 
  shouldn't be affected by survival, but its possible that better growth leads 
  to more mRNA which leads to more accurate quantitation.

Next steps:

- Measure standard curves for the primers I used.

- Repeat the experiment with replicates and more controls, but perhaps only 
  with sz224.

- Repeat the experiment with DNase, and with the −RT control.

- Repeat experiment targeting different regions of the reporter gene.

- Measuring gene expression using the HIS/URA plasmids I've used in previous 
  B1H experiments.  I could even include the two-plasmid system.  I think this 
  will be a valuable point of reference, because these GFP/RFP plasmids are 
  quite different than the HIS/URA ones, and it'd be good to know if the 
  differences have any noticeable impact.

  I'll need to pick primers that target HIS/URA and AmpR (as the house-keeping 
  gene).

- I wonder if using specific primers (instead of random hexamers) for the RT 
  reaction would give better signal.

  - Read about reverse transcriptases
  - Order enzyme that would work well for my purposes.
