#!/usr/bin/env bash
set -euo pipefail

sw warn "Make sure I have 3-12% native PAGE gels in stock." |

# repA gene:
# - Should I include mWasabi (f66) or not (f63)?  Probably I should; I want to 
#   know where the protein is after all.
#
# - I thought about using f17 instead of f66.  It's the same gene but without 
#   the buffer, so it's 7% shorter.  Not a big difference, but potentially 
#   helpful at the margins.  I think f66 might has better primers, though, so 
#   I'm going to stick with it.
#
# Improve PCR yield:
# - NEB suggests ways to improve PCR yield:
#   https://www.neb.com/faqs/2012/01/02/why-do-i-have-a-low-yield-of-pcr-product1
#
#   - increase primer concentration from 0.2 µM to >0.5 µM.
#   - increase cycles to 45
#
# - Note that NEB's Q5 protocol already calls for 0.5 µM primers.
#
# PCR volume:
# - From the :doc:`appendix`, the concentration of ribosomes in a PURExpress 
#   reaction is 2 µM, and it is estimated that each ribosome is recycled 5 
#   times.
#
# - I want to exceed the ribosome concentration, so let's conservatively aim 
#   for 4 µM.
# 
# - I can add at most 1.4 µL DNA to a 5 µL PURExpress reaction, so I would need 
#   14.3 µM stock DNA to reach a final concentration of 4 µM.
# 
# - The MW of f66 is 1441.4 kDa, so 14.3 µM corresponds to 20.6 µg/µL.
#
# - Assuming that I get 100 ng/µL from a PCR reaction (see :expt:`129`), I 
#   would need to concentrate such a reaction 206x.
# 
# - I'll need 1.4 µL for each PURExpress reaction.  I'll need more than that to 
#   do a serial dilution; probably at least 4 µL.  
#
# - That corresponds to an 824 µL PCR reaction.
# 
# - My Qiagen spin cleanup columns have a capacity of 5 µg.  I'm trying to make 
#   82.4 µg of DNA, so I'd need to use ≈16 spin columns for this reaction.
# 
# - That's a lot.  I could do it, though.
#
# - It might be better to come up with some scheme based on plasmid DNA, e.g.:
#
#   - Label via nicking.
#   - Mix miniprep/digest DNA with PCR DNA.
#     - So not all DNA is labeled, but that's fine.
# 
# Negative control:
# - Not sure if I need one.
# - I know repA works; I'm not trying to test that.
# - I'm looking for DNA-concentration-dependent changes in mWasabi expression 
#   and/or migration.
# - My control could just be a lane with ≈75 nM template stock conc.
sw pcr -u f66 -v 500 |

sw spin_cleanup |

sw lyo -v 4 µL |

# End up near 75 nM...
sw serial 4 µL 4µM / 2 7 -0 |

sw ivtt |

sw gel


