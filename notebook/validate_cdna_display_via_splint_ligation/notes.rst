*****************************************
Validate cDNA display via splint-ligation
*****************************************

I'm suspicious that RNase H contamination in PURExpress is preventing 
:expt:`53` from working, due to the RNA-DNA duplex inherent to Y-ligation.  One 
possible way to get around this is to replace Y-ligation with splint ligation.

Splint ligation is the technique normally used in mRNA display to attach 
puromycin to the mRNA.  [Naimudden2016]_ used Y-ligation instead for cDNA 
display, claiming that the latter is more efficient.  In order to use splint 
ligation for cDNA display, I'd need to modify the oligos to add an RT primer, 
but this should be doable.

My thought is to begin by trying to replicate [Barendt2013]_, and then to worry 
about adding the RT primer.
