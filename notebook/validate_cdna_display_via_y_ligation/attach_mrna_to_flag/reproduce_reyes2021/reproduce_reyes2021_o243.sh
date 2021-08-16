#!/usr/bin/env bash
set -euo pipefail

sw zap |

# Controls:
# - Based on expt 101: Without SDS PAGE
sw cond reproduce_reyes2021_cond.xlsx |

# f127 concentration:
# - 363 nM; same as in IVTT reaction
# - See below for discussion.
#
# f127 volume:
# - Total volume of 5 µL to load 2.5 µL of each control.
# - I only really need 4.2 µL, because the +salt control will have some volume 
#   of salt added.  But it's not bad to have some extra.
sw step "Setup the −PUREfrex control:
    ~0.91 µL 2 µM f127
    ~4.09 µL nuclease-free water
    ~Include this control in the 37°C incubation." |

# IVTT kit:
# - I'm still not sure whether PUREfrex 1.0 or 2.0 is better
# - I'm using 1.0 here because I'm trying to reproduce [Reyes2021].
#
# Reaction volume:
# - Need at least 2.5 µL of material (after adding salt) to run the gel.
# - To help make pipetting more accurate, I increased the volume as much as I 
#   could without using more than one aliquot.
#   - 1 aliquot: 12.1 µL (see `sw aliquot_purefrex1`)
#   - I want to use 11 µL for my master mix, so that I have a pipetting buffer 
#     of 10% of the aliquot.
#   - 11 µL / 6.6x × 2 = 3.33 µL reactions
# - If I add more samples in the future, I might have to decrease the reaction 
#   volume accordingly.
#
# mRNA concentration:
# - [Reyes2021]:
#   - 2 µg mRNA in 100 µL reaction
#   - Using the MW of f112 (55110 Da), this is: 363 nM
#
# - expt 99:
#   - 1 µM is ideal (in my hands)
#   - With 2 µM stock, I can only reach 660 nM.
#
# - In the interest of reproducing [Reyes2021], I'm going to use 363 nM.
#
# Reaction time:
# - [Reyes2021] calls for 30 min.
#
# FLAG visualization:
# - o243 is labeled with FITC, so I won't be able to distinguish the RNA 
#   species from any peptide labeled with FluoroTect.
# - In my previous experiment with o236, I didn't observe any FLAG expression 
#   in the TBE/urea gel using FluoroTect.  Maybe there just wasn't any 
#   expression, or maybe the urea gel just didn't resolve it from the 
#   unincorporated FluoroTect.
# - [Reyes2021] only observed the FLAG peptide by Western blot.  In most 
#   experiments, they measured coupling only on on the basis of the difference 
#   between the conjugated and unconjugated mRNA bands (Fig 2).
# - Based on all this, I'm going to not use FluoroTect in this reaction.
sw ivtt \
  f127 f111 \
  -p frex1/flag \
  -d 'FluoroTect GreenLys' \
  -n 6 \
  -v 3.333 \
  -C 2000 \
  -c 363 \
  -r \
  -t '30 min' |

# Salt concentrations:
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

sw step "Setup 8 coupling reactions:
    ~3.33 µL translation reaction
    ~1.53 µL salt dilution" |

sw step "Setup the +salt control:
    ~2.00 µL −PUREfrex control
    ~0.92 µL 3.17x salt solution" |

# Incubation time and temperature from [Reyes2021].
sw step "Incubate at 37°C for 1h." |

# Gel chemistry:
# - [Reyes2021] used SDS/urea PAGE, which I don't have.
#
# - TBE/urea PAGE:
#   - Don't need desalting.
#     - But might want it anyways.
#   - FLAG has pI=4.85, so it doesn't need SDS to run in the right direction.
#   - FLAG also doesn't need a reducing agent, because it doesn't have any 
#     cysteines.
#   - Would probably resolve the mRNA better, and the mRNA is what will give me 
#     the best view of the coupling:
#     - Cy5 is bright.
#     - I'll be able to see both coupled and uncoupled mRNA.
#     - The FLAG might be hard to see; more interference.
#   - Might better denature the ribosomes.
#   - Won't be able to use ladder.
#
# - SDS PAGE:
#  - Would probably resolve the FLAG peptide better, although:
#    - I didn't test TBE/urea in expt 117.
#    - The FLAG peptide will be retained by the desalting column anyways.
#  - Would denature the IVTT proteins, but not sure if that would have any 
#    effect on the gel.
#
# - I think TBE/urea PAGE is the best alternative.
sw gel urea/ivtt 8 -S | sw laser blue
