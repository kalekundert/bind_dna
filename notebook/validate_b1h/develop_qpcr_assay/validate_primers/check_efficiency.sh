#!/usr/bin/env bash
set -euo pipefail

# Template:
# - p221-p228 all have the amplicons
# - f177 is the original gBlock
#
#
# - Serial dilution of plasmid template
#
# - 
# Dilutions:
# - 7x

# https://www.sigmaaldrich.com/US/en/technical-documents/protocol/genomics/qpcr/qpcr-efficiency-determination 
sw serial 200ng/µL / 10 -n 7 -v '20µL' -0 -m 'p221' -d 'water' |
sw pcr -p ssoadv p221,sr1,sr2 -n 8 -d 3 -m primers --skip-thermocycler |
sw pcr -p ssoadv p221,sr3,sr5 -n 8 -d 3 -m primers

