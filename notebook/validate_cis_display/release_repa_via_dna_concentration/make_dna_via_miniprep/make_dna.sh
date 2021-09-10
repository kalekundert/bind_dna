#!/usr/bin/env bash
set -euo pipefail

# Miniprep vs. midiprep
# =====================
# - Plasmids:
#
#   - p124:
#       
#     - trimmed buffer : T7 : RBS : StrepTag : mWasabi : repA : CIS : oriR
#     - MW: 2548181.8
#     - length: 4124
#     - amplicon: f66
#
#   - p27:
#     
#     - T7 : RBS : StrepTag : mWasabi : repA : CIS : oriR
#     - MW: 2378880.9
#     - length: 3850
#     - amplicon: f17
#
# - From :expt:`131`, I want 4 µL of 14.3 µM DNA.
#
# - That works out to 36.4 µg/µL, or 145.8 µg for 4 µL.
#
# - According to Qiagen, I can get "up to" the following yields:
#
#   https://www.qiagen.com/us/products/discovery-and-translational-research/dna-rna-purification/dna-purification/plasmid-dna/qiagen-plasmid-kits/
#
#   - miniprep: 20 µg
#   - midiprep: 100 µg
#
#   These numbers are probably a bit inflated.  I consider 12.5 µg to be an 
#   average miniprep yield (50 µL at 250 ng/µL).
#
# - That works out to:
#
#   - 7-12 minipreps
#   - 2-3 midipreps
#
# - Let me start with two midipreps
#
#   - In the past my midipreps have given almost no yield, but hopefully I'll 
#     have better luck this time.

sw miniprep 

# Plasmid digestion
#
# - p124:
#
#   - AvrII: close to o189
#   - ZraI: close to o186, blunt
#   - AatII: cloe to o186
#
#   Digestion product will be missing T7 terminator...
