#!/usr/bin/env bash
set -euo pipefail

# Bead quantity:
# - Assume I want 25 µL of ligated mRNA stock.
# - IVTT reactions require 1 µM mRNA.
# - That requires a stock concentration of ≈2 µM.
#   - 25 µL @ 2 µM = 50 pmol
# - Only about half the mRNA is ligated.
#   - 100 pmol
# - 1 mg (250 µL, $17) streptavidin beads (NEB S1420) binds >500 pmol 
#   biotinylated 25 nt ssDNA.

# Elution
# - IDT has some good ideas about elution:
#   
#   - DesthioBiotin-TEG: doesn't bind as tightly as regular biotin, so can be 
#     eluted with regular biotin
#
#   - PC Biotin: photocleavable (UV).  Leaves phosphate group.
#
#
