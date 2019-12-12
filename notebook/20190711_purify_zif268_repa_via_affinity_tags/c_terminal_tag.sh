#!/usr/bin/env bash
set -euo pipefail

klab_mutagenesis \
  CTERM_STREP \
  acgtaaccgcaattacagccggctggccacagcttctccctgaaagtgatctcctcagaataatccggcctgcgccggag \
  acgtaaccgcaattacagccggctggccacagcttctcccTGGAGCCACCCGCAGTTCGAAAAAtgaaagtgatctcctcagaataatccggcctgcgccggag \

klab_mutagenesis \
  CTERM_TWIN \
  acgtaaccgcaattacagccggctggccacagcttctcccTGGAGCCACCCGCAGTTCGAAAAAtgaaagtgatctcctcagaataatccggcctgcgccggag \
  acgtaaccgcaattacagccggctggccacagcttctcccAGCGCATGGAGTCATCCTCAATTCGAGAAAGGTGGAGGTTCTGGCGGTGGATCGGGAGGTTCAGCGTGGAGCCACCCGCAGTTCGAAAAAtgaaagtgatctcctcagaataatccggcctgcgccggag \

klab_mutagenesis \
  NTERM_STREP \
  tacgactcactatagggcttaagtataaggaggaaaaaatatgagagacctacaagacgcggtctcagggggaggaggatcag \
  tacgactcactatagggcttaagtataaggaggaaaaaatatggctagcgcaTGGAGTCATCCTCAATTCGAAAAAGGCGCCagagacctacaagacgcggtctcagggggaggaggatcag \

klab_mutagenesis \
  NTERM_TWIN \
  tacgactcactatagggcttaagtataaggaggaaaaaatatggctagcgcaTGGAGTCATCCTCAATTCGAAAAAagagacctacaagacgcggtctcagggggaggaggatcag \
  tacgactcactatagggcttaagtataaggaggaaaaaatatggctAGCGCATGGAGTCATCCTCAATTCGAGAAAGGTGGAGGTTCTGGCGGTGGATCGGGAGGTTCAGCGTGGAGCCACCCGCAGTTCGAAAAAGGCGCCagagacctacaagacgcggtctcagggggaggaggatcag \



klab_mutagenesis \
  CTERM_HIS \
  acgtaaccgcaattacagccggctggccacagcttctccctgaaagtgatctcctcagaataatccggcctgcgccggag \
  acgtaaccgcaattacagccggctggccacagcttctcccGAAAACCTGTACTTCCAATCCCACCACCATCACCATCATtgaaagtgatctcctcagaataatccggcctgcgccggag \

klab_mutagenesis \
  NTERM_HIS \
  tacgactcactatagggcttaagtataaggaggaaaaaatatgagagacctacaagacgcggtctcagggggaggaggatcag \
  tacgactcactatagggcttaagtataaggaggaaaaaatatgCACCACCATCACCATCATGAAAACCTGTACTTCCAATCCGGCGCCagagacctacaagacgcggtctcagggggaggaggatcag \


