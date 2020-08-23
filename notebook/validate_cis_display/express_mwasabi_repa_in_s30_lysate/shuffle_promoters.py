#!/usr/bin/env python3

import random
random.seed(0)

def shuffle_seq(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return ''.join(seq_list)

shuffled_tac = shuffle_seq('TTGACAATTAATCATCGGCTCGTATAATG')
shuffled_t7 = shuffle_seq('TAATACGACTCACTATAGGG')

if __name__ == '__main__':
    print("Shuffled t7:  ", shuffled_t7)
    print("Shuffled tac: ", shuffled_tac)
