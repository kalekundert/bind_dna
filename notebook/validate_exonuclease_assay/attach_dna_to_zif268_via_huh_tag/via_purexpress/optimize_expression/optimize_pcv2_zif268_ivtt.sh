#!/usr/bin/env bash
set -euo pipefail

sw zap |

sw cond optimize_pcv2_zif268_ivtt_cond.xlsx |

# DNA concentrations:
#
# - When I've done this experiment previously (:expt:`99`), I was using RNA 
#   templates, so I can't just use the same concentrations here.
#
# - NEB recommends 25-1000 ng for a 25 µL reaction (5-200 ng for a 5 µL 
#   reaction).  Converting to units of molarity (p107 MW=4089393.4 Da) and 
#   rounding slightly, the recommended final conentrations are: 0.25-10 nM
#
# - Previously I tested 100x concentration ranges in 7 steps (:expt:`99`), and 
#   that worked well.  The recommended range is narrower than that, so I'll 
#   have to expand it.  I want to be sure to saturate the reaction, so I'll 
#   expand more on the larger end: 0.2-20 nM
#
# - I'd like to add 1 µL of DNA to each 5 µL reaction.  That means I need to 
#   reach a stock template concentration of 5×20=100 nM.  If necessary I can 
#   use a slightly lower concentration by adding up to 1.4 µL template to each 
#   reaction.

sw serial 2.6 100nM to 1 7 -0 -m p107 -d 'nuclease-free water' |

sw ivtt p107 -v 5 -n 8 -c 20 -C 100 |

sw gel bolt/ivtt/pcv-zif 8

