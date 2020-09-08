#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw custom "Allow about 6-7h for the entire protocol." |
sw zap |
sw samples.txt |
sw cdna/anneal 1 f85 o129 -v 6 -R 5 | 
sw custom "Dilute the following mRNA samples to 333 nM using nuclease-free water:~5 µM mRNA (f85): dilute 1 µL to 15 µL~1.25 µM annealed but unligated mRNA/linker: dilute 1 µL to 3.75 µL" |
sw cdna/ligate 1 -v 42.5 |
sw custom "Reserve 2.5 µL of the ligated but unfiltered mRNA/linker." |
sw cdna/wash |
sw purex "mRNA/linker" "−ligation" "−linker" -T 333 -r -v 3.6 |
sw note $'Target mRNA concentration: #18\nStock mRNA concentration: #19\nReaction volume: Require only 1 PURExpress aliquot (5 µL)' "reactions" -W |
sw gel sds/o194 o129,a,al,alf,alfx,ax,x−o129,f85 -S |
sw laser blue red
