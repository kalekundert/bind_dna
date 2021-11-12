#!/usr/bin/env bash
set -euo pipefail

# Estimated time required: 4h30

sw zap |

sw cond attach_via_purefrex_cond.xlsx |

# f89 volume:
# - Total volume of 5 µL to load 2.5 µL of each control.
# - I only really need 4.2 µL, because the +salt control will have some volume 
#   of salt added.  But it's not bad to have some extra.
sw step "Setup the −PUREfrex control:~1.90 µL 1 µM f89~3.10 µL nuclease-free water~Include this control in the 37°C incubation." |

# reaction time:
# - Most protocols have 30 min incubations.
# - [Reyes2021] shows that the reaction is basically complete at 5 min.
# - But I'm making mWasabi, and most of the other papers are making peptides.  
#   I've also seen pretty low levels of mWasabi fluorescence after 30 min, so I 
#   think it's reasonable to give it more time
#
# reaction volume:
# - Need at least 2.5 µL of material (after adding salt) to run the gel.
# - To help make pipetting more accurate, I increased the volume as much as I 
#   could without using more than one aliquot (see `sw aliquot_purefrex1`).
# - If I add more samples in the future, I might have to decrease the reaction 
#   volume accordingly.
sw ivtt f85 f89 -p frex1/gfp -n 5 -v 4 -x 20 -r -t 0 |
sw ivtt f85 f89 -p frex2/gfp -n 5 -v 4 -x 20 -r -t '60 min' |

# salt concentrations:
# - Want to reproduce the titration in [Reyes2021] (Fig 2a).
#
#   - The reaction seems more sensitive to the salt concentration than to any 
#     other parameter, and the optimal concentration might be different for 
#     each mRNA.  I'm hoping that a titration will give me the best chances of 
#     seeing some successful coupling.
#
#   - The highest concentration in the [Reyes2021] titration (750 mM KCl, 65 mM 
#     MgOAc) is actually taken from [Naimudden2016], which is incidentally the 
#     concentration that I've used previously.
#
# - Given these target concentrations, I used the `cdna/couple` protocol to 
#   calulate the corresponding volumes:
#
#     $ sw cdna/couple -n 2 -v 4 -k 750 -m 65
sw step "Prepare a 3.17x salt solution:~5.00 µL 1M MgOAc~19.23 µL 3M KCl~Keep at room temperature." |

# Same dilution as [Reyes2021] (Fig 2a), except without the 1/8 dilution 
# because I only want to use 1 aliquot of PUREfrex.
sw serial 10µL 3.17x / 2 3 -0 -m "salt solution" -d "nuclease-free water" |

# Volume of salt dilution from above `sw cdna/couple` command.
sw step "Setup 8 coupling reactions:~4.00 µL translation reaction~1.84 µL salt dilution" |

sw step "Setup the +salt control:~2.0 µL −PUREfrex control~0.92 µL 3.17x salt solution" |

sw step "Incubate at 25°C for 1h." |

sw gel urea/ivtt 12

