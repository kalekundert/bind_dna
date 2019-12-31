***************************
Optimize mRNA transcription
***************************

The mRNA in the 2019/12/9 experiment of 
:expt:`20191209_ligate_linker_n_via_naimudden2016` looked noticeably degraded.  
There are a few reasons why this may have happened:

- Too much RNase (leftover from miniprep) in template DNA.  I could test for 
  this by doing a side-by-side comparison between the template I used this time 
  and template that I purified by phenol-chloroform extraction.

- Expired HiScribe kit.  I remember getting smeary RNA with older IVT kits 
  in grad school, so it wouldn't surprise me if that was happening here.  If 
  none of my other ideas explain the heterogeneity, I should just try 
  ordering a new kit.

- I contaminated something.  I could test for this by repeating the 
  transcription and just being more careful.  I thought I was pretty careful 
  this time, though...

- I left the mRNA at 4°C overnight, because I started the in vitro 
  transcription reaction right before going home.  

I want to see I can understand what happened and improve the quality of my 
mRNA.

Results
=======

PCR & phenol/chloroform --- 2019/12/16
--------------------------------------
In this experiment, I want to test if I can reduce degradation by better 
purifying the template DNA.  If this works, it would imply that the degradation 
in the 2019/12/9 experiment is due to RNase activity in the template (which is 
plausible, since the template was basically miniprepped DNA).

I'm comparing two ways of purifying the template:

- PCR: Basically this just lets me really dilute the plasmid template, which 
  may be contaiminated with RNase A from the miniprep. Note that I still need 
  to do the XmnI digestion, because PCR primers in the Y-tag don't amplify 
  well: :expt:`20191004_linearize_cdna_display_gene`

- Phenol/chloroform extraction: Gold-standard method for removing protein 
  contamination of any kind from nucleic acid samples.

In either case, I purified and concentrated the DNA using magnetic beads before 
in vitro transcription.

.. protocol:: 20191216_pcr.txt 20191216_phenol_chloroform_extraction.txt 
   20191216_transcribe_rna.txt
   
   See binder.

.. figure:: 20191217_compare_template_dna.svg

   In vitro transcribed RNA, with the DNA template purified as indicated.

- All four transcription reactions worked well.  The products also the right 
  size, which I couldn't confirm previously because I didn't have a ladder.

- The −/− template is the same template I used in the 2019/12/9 experiment.  It 
  looks quite clean.  This suggests to me that the degradation in the 2019/12/9 
  experiment was due to neither RNase contamination in the template nor an 
  expired T7 kit.  My hypothesis now is that it was due to letting the reaction 
  sit overnight at 4°C.

Discussion
==========
- No purification beyond using magnetic beads after XmnI digestion seems 
  necessary to get clean transcription.
