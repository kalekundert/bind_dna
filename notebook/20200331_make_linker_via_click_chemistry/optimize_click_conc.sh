#!/usr/bin/env zsh

stepwise direct_dilution 2µL 400µM 200µM 4 -m o127 | stepwise <<EOF
- Run 4 click reactions:

  Reagent       Stock  Volume    4.4x
  ───────────────────────────────────
  o126         400 µM  1.0 µL  4.4 µL
  o127     400–200 µM  2.0 µL
  PBS              4x  1.0 µL  4.4 µL
  ───────────────────────────────────
                               2.0 µL/rxn

  - Incubate for the optimized length of time.

- Run an E-gel:

   - Dilute 0.2 µL of each reaction to 40 µL.
   - Load 20 µL into a 2% E-gel [1].
   - Run for 10 min.

Footnotes:

[1] ≈100 ng/oligo in each lane.
EOF

# vim: tw=53
