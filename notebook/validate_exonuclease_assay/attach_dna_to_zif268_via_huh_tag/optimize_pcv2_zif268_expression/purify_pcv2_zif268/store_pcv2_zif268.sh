#!/usr/bin/env bash
set -euo pipefail

# Storage buffer:
# - Binding buffer with the following modifications:
#   - Without nonspecific binding agents.
#   - With 50% glycerol.
#
# - The following recipe most resembles [Zykovich2009].  I like the recipe from 
#   [Zykovich2009] because all the counter-ions are the same, and there seem to 
#   be no frills.  I replaced KCl with NaCl to make running SDS PAGE gels 
#   easier:
#
#   - 10 mM tris
#   - 90 mM NaCl
#   - 1 mM MgCl₂
#   - 90 µM ZnCl₂
#   - 5 mM DTT
#   - pH=7.5
#
# - I don't have ZnCl₂; I'll just use ZnOAc instead since that's what I've used 
#   before.
#
# - 50 mM salt (or less) is recommended for native PAGE:
#   https://assets.thermofisher.com/TFS-Assets/LSG/manuals/nativepage_man.pdf
#
#   - But I'd rather not adjust the binding buffer to accomodate the EMSA 
#     assay, and I can dilute the smaples before EMSA anyways.
#
#   - Maybe I should reduce [NaCl] to 50 mM?  It's not unprecedented; the 
#     Zif268 binding buffer from [Liu2005] has 50 mM NaCl.
#
# - I'll add the glycerol directly to the 1x buffer above.  Technically that 
#   will make the buffer 0.5x.
#
# - Store at -20°C
#   - It may not work to simply add glycerol and store at -20°C, but I'm just 
#     going to try it and worry about it if I have problems.

# Protein concentration:
#
# - 
#   https://www.sheffield.ac.uk/polopoly_fs/1.779623!/file/PurificationGuide.pdf#
#
#   - "Please note that during centrifugation protein concentration becomes 
#     very high at the bottom of the concentrator which could cause protein 
#     precipitation.  For some proteins it is useful to do concentration in 
#     short (5-10 min) shots with the mixing of the sample after each shot."
#
#   - "For concentrated protein solutions (more than 10 mg/ml) it is not 
#     necessary to add anti-freezing agent such as glycerol, but for dilute 
#     solutions it is better to add glycerol to a 10-40% concentration."
#
# - 
#   http://tools.thermofisher.com/content/sfs/brochures/TR0043-Protein-storage.pdf
#
#   - "Dilute protein solutions (<1 mg/ml) are more prone to inactivation and 
#     loss as a result of low-level binding to the storage vessel.  Therefore, 
#     it is common practice to add “carrier” or “filler” protein, such as 
#     purified bovine serum albumin (BSA) to 1-5 mg/ml (0.1-0.5%), to dilute 
#     protein solutions to protect against such degradation and loss."
#
# - I can also think about what molarity I want in my reactions, and make 10x 
#   aliquots.
#
#   - EMSA:
#     - My previous EMSA experiments used PURExpress, so unknown protein 
#       concentration.
#
#     - My DNA concentration in those experiments was 625 pM, though.  I don't 
#       know the optimal protein:DNA ratio, but 1:1 is a reasonable assumption 
#       (although it is also a lower limit):
#
#       - PCV2-Zif268: 27.9 kDa
#       - 625 pM = 17.4 pg/µL
#       - This is much too little to visualize (I think the detection limit for 
#         most protein stains is 1-10 ng).  But I don't need to visualize the 
#         protein for EMSA.
#
#     - The EMSA I'm planning (`../confirm_zif268_activity`) will require a 
#       stock concentration of at least 1.6 µM:
#
#       - 1.6 µM = 44.6 µg/mL
#
#     - Cannot find recommendation from Thermo/Life Technologies rearding how 
#       much material to load on native PAGE gel (if I want to visualize the 
#       protein).  Bolt gels recommend 250 ng/band, but in my experience this 
#       is much too low.  
#
#   - Exonuclease
#
#     - These reactions can also be very dilute, because I'm detecting with 
#       qPCR and I want to avoid trans-interactions.
#
# - I want my stock to be in units of molarity:
#
#   - 1 mg/mL = 35.8 µM seems like a reasonable/practical target.
#   - Round numbers:
#     - 32 µM = 892.8 µg/mL
#     - 16 µM = 446.4 µg/mL (if I can't reach 1 mg/mL)
#
# I think 1 mg/mL is a reasonable target.  We'll see what I can actually get.

# Amicon filter:
# - Pore size:
#   - PCV2-Zif268: 27.9 kDa
#   - Similar to α-chymotrpsinogen (25 kDa)
#   - https://www.emdmillipore.com/US/en/product/Amicon-Ultra-0.5mL-Centrifugal-Filters-for-DNA-and-Protein-Purification-and-Concentration,MM_NF-C82301#overview
#   - 10K filter is best
#
# - Volume:
#   - Only need 2 mL
#   - But if I decide to buy something, I think 4 mL would be more 
#     future-proof.
#
# - Isaac has a stash of 4 mL, 10K filters.  As long as I don't need too many, 
#   I can use those.

# Buffer exchange:
# - Elution buffer:
#   - 50 mM sodium phosphate
#   - 300 mM NaCl
#   - 0.05% tween
#   - 250 mM imidazole
#   - pH=8.0
#
# - Each spin:
#   - Concentrate sample to 50 µL
#   - Add 4 mL buffer
#   - 80x dilution factor
#
# - After first dilution:
#   - 1 mM buffer --- negligible rel. to storage buffer
#   - 4 mM NaCl --- negligible rel. to storage buffer
#   - .0006% tween --- doesn't seem like a lot, but don't really have a point 
#     of comparison
#   - 3 mM imidazole --- less than wash buffer
#
# - I think 1 dilution will be enough.

sw step "Measure concentration by nanodrop~Extinction coefficient: 17460 M⁻¹cm⁻¹ [ExPASy ProtParam]" |
sw storage_buffer.txt |
sw step "Concentrate the purified protein and tranfser it to storage buffer:
      ~ Load the eluted protein on a 4 mL, 10K Amicon spin filter.
      ~ Spin 7000g, 4°C, 15 min.
      ~ Add 4 mL Zif268 storage buffer.
      ~ Spin 7000g, 4°C, 15 min.
      ~ Dilute to 64 µM (1.785 mg/mL)." |
sw step "Add 1 volume 80% glycerol.  Mix well." |
sw step "Store at −20°C."


