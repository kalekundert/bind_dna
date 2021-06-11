#!/usr/bin/env bash
set -euo pipefail

# Templates:
# - Want to measure yield for both Zif268 and FLAG
#
# Sample quantity:
# - [Nilson2013]: 10-20 µg per lane (They call for an entire 50 µL IVT 
#   reaction.  I looked up the referenced IVT protocol, and this is the 
#   expected yield.)
#
# - My samples are 10 µM, which corresponds to:
#   - f11: 1260 ng/µL
#   - f111: 398 ng/µL
#   
#   I'm loading 5 µL, which corresponds to:
#   - f11: 6 µg
#   - f111: 2 µg
#
#   So I'm a little short of [Nilson2013], but in the right ballpark.
#
# Prerunning:
# - Denaturing gels are pre-run because the tempature helps with denaturing.
# - But in my experience I get perfectly good gels without it.

sw gel urea/purify f11,f111 -S |

sw step "Cut the desired band out of the gel.
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

sw step "Incubate overnight at room temperature with agitation.
~4°C would also be ok.
~thermomixer at 800 rpm or nutator." |

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

sw step "Measure yield by NanoDrop"
