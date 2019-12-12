#!/usr/bin/env bash
set -euo pipefail

klab_mutagenesis \
  TAC \
  CACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACGcgtaatacgactcactatagggCTTAAGTATAAGGAGGAAAAAAtATGGCTAGC \
  CACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACGcgTTGACAATTAATCATCGGCTCGTATAATGCTTAAGTATAAGGAGGAAAAAAtATGGCTAGC \

klab_mutagenesis \
  RBS \
  aagaaacgcgtaatacgactcactatagggCTTAAGTATAAGGAGGAAAAAAtATGgctagctggagccacccgcagttcgaaaaa \
  aagaaacgcgtaatacgactcactatagggAATAGAAAAAtGATCAGTTAAAGGGTgctagctggagccacccgcagttcgaaaaa \

