#!/usr/bin/env bash
set -euo pipefail

# f112, f118: spacer-18 linkers
# f113, f119: poly-A linkers

sw zap |

# mRNA volume:
# - Let's assume I want to end up with 20 µL ligated mRNA stock.
# - IVTT reactions require 1 µM mRNA.
# - Reaching that requires a stock concentration of ≈2 µM.
#   - 20 µL @ 2 µM = 40 pmol
# - ≈40% of the mRNA will be ligated:
#   - 40 pmol / 40% = 100 pmol
# - ≈70% of the mRNA will be recovered from the gel:
#   - 100 pmol / 70% = 142 pmol = 14.2 µL @ 10 µM
#   - f11: 18 µg
#   - f111: 6 µg
# - The recommended amount of material per lane when doing PAGE gel 
#   purification is 10-20 µg, so this is about right.
# - I might even lyophilize to get the volume down to 5 µL, so I can fit the 
#   entire sample in one lane.

# Ligation protocol:
# - I don't think the [Reyes2021]_ protocol is optimal, because it dilutes to 
#   oligos before the annealing reaction.
# - But it's the best protocol I have right now, and I don't expect 
#   optimization to significantly improve yield.
# - In the interest of conserving expensive material, I'm going to use only 
#   0.6x linker.

# mRNA/linker concentrations:
# - [Reyes2021] has final concentrations:
#   - mRNA: 4 µM
#   - linker: 6 µM
# - My stocks are not concentrated enough to reach this, but I can reach:
#   - mRNA: 4 µM
#   - linker 2.4 µM (0.6x)

sw reaction \
  -s "ligation" \
  'nuclease-free water  ;         ; to 4.8 µL ; +' \
  'T4 RNA ligase buffer ;     10x ;    0.5 µL ; +' \
  'ATP                  ;   10 mM ;    0.5 µL ; +' \
  'f111,f11             ;   10 µM ;    2.0 µL ; -' \
  'o236,o128            ;   10 µM ;    1.2 µL ; -' \
  -V 10 \
  -n 2 |

sw step "Incubate as follows:~90°C for 30s~Cool to 25°C at 1°C/s" |

# I'd rather include the enzyme in a master mix (e.g. anneal with just the 
# oligos and 1x T4 RNA ligase buffer; add water, ATP, more T4 RNA ligase 
# buffer, and enzyme after incubation).  This would avoid having to pipet such 
# small volumes, but for now I'm just going to follow [Reyes2021] exactly.
sw step "Add to each reaction:~2 µL 10 U/µL T4 RNA ligase" |

sw step "Incubate at 25°C for 30 min." |

sw step "Label the products: f112, f118" |

sw step "Concentrate the reaction to ≈5 µL by lyophilization." |

sw gel urea/purify f112,f118 -S |

sw step "Cut the ligated and unligated mRNA bands out of the gel.
~Place the gel over a TLC plate.
~Use a UV light to visualize the RNA (dark spot)
~Consider visualizing remaining gel to ensure that all desired RNA was excised." |

sw note "Based on Fitzy's DNA PAGE purification protocol, [Nilson2013], and [Petrov2013]." |

sw step "Crush gel slices.
~Poke 3-4 holes in 0.65 mL tube with 27 g needle.
~Place gel slice inside 0.65 mL tube.
~Place 0.65 mL tube inside 1.5 mL tube.
~Centrifuge 5 min at max speed." |

sw step "Resuspend gel in 400 µL PAGE elution buffer
~Tris-HCl, pH 7.5 10 mM
~NaCl 500 mM
~EDTA 1 mM
~SDS 0.1%" |

sw step "Incubate at 4°C overnight with end-over-end mixing." |

# Filter material:
# - Corning has a good guide on which material to select:
#   https://www.corning.com/catalog/cls/documents/selection-guides/t_filterselectionguide.pdf
#
# - I want 0.22 µm, because that's the standard size for filter sterilizing 
#   biological buffers (that's not my application here, but I can see myself 
#   wanting to do that).
#
# - I want cellulose acetate filters.  Nylon and cellulose nitrate have high 
#   DNA binding, which will cause me to lose material.  The downside to 
#   cellulose acetate is that it has a wetting agent that will end up in the 
#   sample.  However, this will be removed by the Zymo column in the subsequent 
#   step.
#
# - Product number: 8161 (non sterile)
#
# Centrifugation speed:
# - Fitzy's DNA PAGE purification protocol calls for 4 min at 7000 rpm
# - The Corning guide (see above) includes an agarose gel purification 
#   protocol, which calls for 13,000g for 5-20 min.  But this protocol has no 
#   incubation step, so I gather that the spin is supposed to pull the solvent 
#   out of the gel.  I probably don't need to go so fast.
# - But why not go as fast as possible?
sw step "Spin-elute through a Spin-X column or similar to remove gel slices
~4 min at 7000 rpm" |

sw spin_cleanup zymo/rna-clean-conc/5 -s 400 |

sw aliquot 5µL 2µM
