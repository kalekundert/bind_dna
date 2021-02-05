#!/usr/bin/env python3

import random
random.seed(0)

def shuffle_seq(seq):
    return ''.join(random.sample(seq, len(seq)))

def shuffled_gibson_insert(seq_5, seq_to_shuffle, seq_3):
    seq_shuffled = shuffle_seq(seq_to_shuffle)
    return seq_5.lower() + seq_shuffled.upper() + seq_3.lower()

## Shuffle CIS

# 20-40 bp overlaps seem to be a common recommendation, and I get the 
# impression that 20 is a minimum and 40 is more of an ideal length.

cis_insert = shuffled_gibson_insert(
        'taaccgcaattacagccggctggccacagcttctccctga',
        'aagtgatctcctcagaataatccggcctgcgccggaggcatccgcacgcctgaagcccgccggtgcacaaaaaaacagcgtcgcatgcaaaaaacaatctcatcatccaccttctggagcatccgattccccctgtttttaatacaaaatacgcctcagcgacggggaattt',
        'tgcttatccacatttaactgcaagggacttccccataagg',
)
print(f"shuffled_cis_insert\t{cis_insert}")

## Shuffle ORI

# Had to reduce the 3' overlap to 20 bp to avoid a repeated sequence that 
# would've interfered with synthesis.

ori_insert = shuffled_gibson_insert(
        'ctgtttttaatacaaaatacgcctcagcgacggggaattt',
        'tgcttatccacatttaactgcaagggacttccccataaggttacaaccgttcatgtcataaagcgccagccgccagtcttacagggtgcaatgtatcttttaaacacctgtttatatctcctttaaactacttaattacattcatttaaaaagaaaacctattcactgcctgtcctgtggacagacag',
        'gctagcaaaaggccagcaaa',
)
print(f"shuffled_ori_insert\t{ori_insert}")

# Shuffle RBS

rbs_seq = shuffle_seq('CTTAAGTATAAGGAGGAAAAAAtATG')
print(f"shuffled_rbs_seq\t{rbs_seq}")
