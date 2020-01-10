#!/usr/bin/env sh

# template: 27

klab_mutagenesis \
        SHUFFLE_TAC \
        cacatttccccgaaaagtgccacctgacgtctaagaaacgcgTAATACGACTCACTATAGGGcttaagtataaggaggaaaaaatatggctagc \
        cacatttccccgaaaagtgccacctgacgtctaagaaacgcgAAGTATTAATTACGCCATGAGCGTTATCTcttaagtataaggaggaaaaaatatggctagc \

klab_mutagenesis \
        SHUFFLE_T7 \
        cacatttccccgaaaagtgccacctgacgtctaagaaacgcgTAATACGACTCACTATAGGGcttaagtataaggaggaaaaaatatggctagc \
        cacatttccccgaaaagtgccacctgacgtctaagaaacgcgAGAATTCAGCTACGATGCTActtaagtataaggaggaaaaaatatggctagc \


