#!/usr/bin/env zsh
set -euo pipefail
alias sw=stepwise

sw zap |
sw cdna/anneal 1 f85 o129 -v 20 -R 5 |
sw custom "Reserve 0.44 µL of the annealed product for the '−ligate' lane of the gel, and dilute it to 4.44 µL with nuclease-free water." |
sw cdna/ligate 1 -v 180 |
sw custom "Reserve 4.44 µL of the ligated product for the '−filter' lane of the gel." |
sw echo "
- Split the ligation reaction into 4 40 µL samples, 
  to be prepared as follows:
  
  urea:  −  +  −  +
  o194:  −  −  +  +

- Remove unligated linker using a MWCO spin filter:
  
  - For samples with o194:
    - Add 5 µL 100 µM o194.
    - Incubate at 95°C for 2 min [1]

  - Bring reactions to 500 µL with 8M urea or water.
  - Load onto a 100 kDa MWCO spin-filter [2].
  - Spin 14000g, 15 min, 4°C.
  - Wash with 500 µL 8M urea or water.
  - Wash with 500 µL water.

  - For samples with urea:
    - Wash with water again.
    - Wash with water again.

  - Invert the filter into a clean tube and spin 
    1000g, 2 min to collect ligated product in a 
    volume of ≈15 µL.

Notes:
[1] When preparing samples for PAGE, this incubation 
    step is 70°C for 3 min.  I chose to use a higher 
    temperature here because there is no formamide in 
    the buffer to help lower melting temperatures.

[2] Amicon UFC510024
" |
sw gel urea/o194 −ligate,−filter,−/−,+/−,−/+,+/+
