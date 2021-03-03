#!/usr/bin/env bash
set -euo pipefail

sw zap |

sw cond attach_via_purefrex_cond.xlsx |

# reaction time:
# - Most protocols have 30 min incubations.
# - [Reyes2021] shows that the reaction is basically complete at 5 min.
# - But I'm making mWasabi, and most of the other papers are making peptides.  
#   I've also seen pretty low levels of mWasabi fluorescence after 30 min, so I 
#   think it's reasonable to give it more time
sw ivtt f89 -p frex/gfp -t '60 min' -v 8*2.5*1.1 -r |

# f89 volume:
# - Chose 2.5 µL to use remainder of thawed aliquots.
sw step "Setup the −PUREfrex control:~2.50 µL 1 µM f89~5.05 µL nuclease-free water" |

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
#     $ sw cdna/couple -n 2 -v 5 -k 750 -m 65
sw step "Prepare a 3.17x salt solution:~5.00 µL 1M MgOAc~19.23 µL 3M KCl" |

# Same dilution as [Reyes2021] (Fig 2a), except without the 1/8 dilution 
# because I only want to use 1 aliquot of PUREfrex.
sw serial 10µL 3.17x / 2 3 -0 -m "salt solution" -d "nuclease-free water" |

# Volume of salt dilution from above `sw cdna/couple` command.
sw step "Setup 4 coupling reactions:~5.0 µL translation reaction~2.3 µL salt dilution" |

# incubation time:
# - The one example I have of a −20°C incubation [Cotten2011] does both the 
#   25°C step and the −20°C step.  So I'll do the same here.
sw step "Incubate at 25°C for 1h." |

sw gel bolt/gfp-mrna 5 |

sw step "Incubate the remaining coupling reaction overnight at −20°C, and run the same gel as above the next day."

