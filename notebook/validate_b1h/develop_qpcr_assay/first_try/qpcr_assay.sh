#!/usr/bin/env bash
set -euo pipefail

sw extract_rna.sh |
sw reverse_transcribe.sh |

# Master mix:
# - I could either put the template or the primers in the master mix.
# - I think I have to go for the templates.  I want to compare the ratio of GFP 
#   to RFP for each sample, and for this it's essential that both tubes have 
#   the same amount of template.
sw pcr -p ssoadv cDNA,sr1,sr2 cDNA,sr3,sr5 -n9 -d2 -m '' |
sw sub 'qPCR reactions' 'qPCR reactions for each sample'
