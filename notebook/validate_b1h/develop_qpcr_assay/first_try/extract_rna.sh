#!/usr/bin/env bash
set -euo pipefail

# Protocol based on:
# ~sgrna/notebook/20180925_quantify_sgrna_levels/quantify_sgrna_levels.txt

sw step "Grow 1 mL overnight cultures of sz221-sz228 in LB+Carb at 37°C for 16h." |
sw note "The small culture size is meant to reduce the chance of cross-contamination when using 24-well blocks." |

# Growth time:
# - There needs to be enough time to:
#   - Express rpoZ-Zif268 (induced with IPTG).
#   - Allow mRNA levels to reach steady state.
#   - Have enough cells to be able to detect mRNA.
#
# - mRNA has a very short half-life, so I think that a stready state will be 
#   reached very quickly and that the most important consideration will just be 
#   the number of cells.
#
# - I probably still want cells in log phase, just because who knows what's 
#   happening in lag phase.  Early log phase should be a good target.
#
# - Right now I don't know the relationship between OD and growth phase.  In my 
#   sgRNA project, I just grew the cells for 2h.  I think that's reasonable for 
#   now, and I can revisit this once I have real OD numbers.
sw step "Inoculate day cultures:~1 mL LB + Carb + 10 µM IPTG~10 µL overnight culture" |
sw step "Grow day cultures at 37°C for 2h." |

sw step "Extract total cellular RNA:
  ~Pellet cells at 7000g for 3 min
  ~Resuspend each cell pellet in 1 mL TRIzol (Invitrogen 15596026).
  ~Incubate 5 min at room temperature.
  ~Add 200 µL chloroform.
  ~Vortex vigorously.
  ~Centrifuge for 15 min at 20,000g and 4°C.
  ~Transfer aqueous phase (top, not pink, ≈500 µL) for each sample to a clean tube, taking care to avoid transferring any of the organic phase." |

sw step "Concentrate and purify the RNA by ethanol precipitation:
   ~Add 1 µL GlycoBlue (15 mg/mL).
   ~Add 500 µL isopropanol.
   ~Incubate at room temperature for 10 min.
   ~Pellet for 20 min at 12,000g and 4°C.
   ~Carefully pipet off all supernatant.
   ~Resuspend pellet in 70% EtOH.
   ~Vortex briefly
   ~Pellet for 5 min at 7,500g and 4°C.
   ~Carefully pipet off all supernatant.
   ~Air dry for 10 min.
   ~Resuspend RNA in 10 µL water." |

sw step "Measure the RNA concentration of each sample by Nanodrop."
