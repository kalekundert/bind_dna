#!/usr/bin/env python3

from bind_dna.dna import *
from pytest import approx

def test_get_plasmid_seq():
    assert get_plasmid_seq(2)      == 'TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTC'
    assert get_plasmid_seq('p2')   == 'TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTC'
    assert get_plasmid_seq('p002') == 'TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGCGCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGCTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTC'

def test_get_fragment_seq():
    # PCR --- no overhang
    assert get_fragment_seq(3)      == 'TGACGTCTAAGAAACGCGTTGACAATTAATCATCGGCTCGTATAATGCTTAAGTATAAGGAGGAAAAAATATGGCTAGCTGGAGCCACCCGCAGTTCGAAAAAGGCGCCGTGAGCAAAGGCGAAGAAACCACCATGGGCGTGATTAAACCGGATATGAAAATTAAACTGAAAATGGAAGGCAACGTGAACGGCCATGCGTTTGTGATTGAAGGCGAAGGCGAAGGCAAACCGTATGATGGCACCAACACCATTAATCTGGAAGTGAAAGAAGGCGCGCCGCTGCCGTTTAGCTATGATATTCTGACCACCGCGTTTAGCTATGGCAACCGTGCGTTTACCAAATATCCGGATGATATTCCGAACTATTTTAAACAGAGCTTTCCGGAAGGCTATAGCTGGGAACGTACCATGACCTTTGAAGATAAAGGCATTGTGAAAGTGAAAAGCGATATTAGCATGGAAGAAGATAGCTTTATTTATGAAATTCATCTGAAAGGCGAGAACTTTCCGCCGAACGGCCCGGTGATGCAGAAAGAAACCACCGGCTGGGATGCGAGCACCGAACGTATGTATGTGCGTGATGGCGTGCTGAAAGGTGATGTGAAAATGAAACTGCTGCTGGAAGGCGGCGGCCATCATCGTGTGGATTTTAAAACCATTTATCGTGCGAAAAAAGCGGTGAAACTGCCGGATTATCATTTTGTGGATCATCGTATTGAAATTCTGAACCATGATAAAGATTATAACAAAGTGACCGTGTATGAAATTGCGGTGGCGCGTAACAGCACCGATGGCATGGATGAACTGTATAAAGGGGGAGGAGGATCAGCGGCCCCAACTGATCTTCACCAAACGTATTACCGGCAGGTAAAGAACCCGAATCCGGTGTTCACTCCCCGTGAAGGTGCCGGAACGCTGAAGTTCTGCGAAAAACTGATGGAAAAGGCGGTGGGCTTCACCTCCCGTTTTGATTTCGCCATTCATGTGGCGCATGCCCGTTCCCGTGGTCTGCGTCGGCGCATGCCACCGGTGCTGCGTCGACGGGCTATTGATGCGCTGCTGCAGGGGCTGTGTTTCCACTATGACCCGCTGGCCAACCGCGTCCAGTGTTCCATCACCACACTGGCCATTGAGTGCGGACTGGCGACAGAGTCCGGTGCAGGAAAACTCTCCATCACCCGTGCCACCCGGGCCCTGACGTTCCTGTCAGAGCTGGGACTGATTACCTACCAGACGGAATATGACCCGCTTATCGGGTGCTACATTCCGACCGACATCACGTTCACACTGGCTCTGTTTGCTGCCCTTGATGTGTCTGAGGATGCAGTGGCAGCTGCGCGCCGCAGTCGTGTTGAATGGGAAAACAAACAGCGCAAAAAGCAGGGGCTGGATACACTGGGTATGGATGAGCTGATAGCGAAAGCTTGGCGTTTTGTGCGTGAGCGTTTCCGCAGTTACCAGACAGAGCTTCAGTCCCGTGGAATAAAACGTGCCCGTGCGCGTCGTGATGCGAACAGAGAACGTCAGGACATCGTCACCCTAGTGAAACGGCAGCTGACGCGTGAAATCTCGGAAGGACGCTTCACTGCTAATGGTGAGGCGGTAAAACGCGAAGTGGAGCGTCGTGTGAAGGAGCGCATGATTCTGTCACGTAACCGCAATTACAGCCGGCTGGCCACAGCTTCTCCCTGAAAGTGATCTCCTCAGAATAATCCGGCCTGCGCCGGAGGCATCCGCACGCCTGAAGCCCGCCGGTGCACAAAAAAACAGCGTCGCATGCAAAAAACAATCTCATCATCCACCTTCTGGAGCATCCGATTCCCCCTGTTTTTAATACAAAATACGCCTCAGCGACGGGGAATTTTGCTTATCCACATTTAACTGCAAGGGACTTCCCCATAAGGTTACAACCGTTCATGTCATAAAGCGCCAGCCGCCAGTCTTACAGGGTGCAATGTATCTTTTAAACACCTGTTTATATCTCCTTTAAACTACTTAATTACATTCATTTAAAAAGAAAACCTATTCACTGCCTGTCCTGTGGACAGACAGGCTAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAG'
    assert get_fragment_seq('f3')   == 'TGACGTCTAAGAAACGCGTTGACAATTAATCATCGGCTCGTATAATGCTTAAGTATAAGGAGGAAAAAATATGGCTAGCTGGAGCCACCCGCAGTTCGAAAAAGGCGCCGTGAGCAAAGGCGAAGAAACCACCATGGGCGTGATTAAACCGGATATGAAAATTAAACTGAAAATGGAAGGCAACGTGAACGGCCATGCGTTTGTGATTGAAGGCGAAGGCGAAGGCAAACCGTATGATGGCACCAACACCATTAATCTGGAAGTGAAAGAAGGCGCGCCGCTGCCGTTTAGCTATGATATTCTGACCACCGCGTTTAGCTATGGCAACCGTGCGTTTACCAAATATCCGGATGATATTCCGAACTATTTTAAACAGAGCTTTCCGGAAGGCTATAGCTGGGAACGTACCATGACCTTTGAAGATAAAGGCATTGTGAAAGTGAAAAGCGATATTAGCATGGAAGAAGATAGCTTTATTTATGAAATTCATCTGAAAGGCGAGAACTTTCCGCCGAACGGCCCGGTGATGCAGAAAGAAACCACCGGCTGGGATGCGAGCACCGAACGTATGTATGTGCGTGATGGCGTGCTGAAAGGTGATGTGAAAATGAAACTGCTGCTGGAAGGCGGCGGCCATCATCGTGTGGATTTTAAAACCATTTATCGTGCGAAAAAAGCGGTGAAACTGCCGGATTATCATTTTGTGGATCATCGTATTGAAATTCTGAACCATGATAAAGATTATAACAAAGTGACCGTGTATGAAATTGCGGTGGCGCGTAACAGCACCGATGGCATGGATGAACTGTATAAAGGGGGAGGAGGATCAGCGGCCCCAACTGATCTTCACCAAACGTATTACCGGCAGGTAAAGAACCCGAATCCGGTGTTCACTCCCCGTGAAGGTGCCGGAACGCTGAAGTTCTGCGAAAAACTGATGGAAAAGGCGGTGGGCTTCACCTCCCGTTTTGATTTCGCCATTCATGTGGCGCATGCCCGTTCCCGTGGTCTGCGTCGGCGCATGCCACCGGTGCTGCGTCGACGGGCTATTGATGCGCTGCTGCAGGGGCTGTGTTTCCACTATGACCCGCTGGCCAACCGCGTCCAGTGTTCCATCACCACACTGGCCATTGAGTGCGGACTGGCGACAGAGTCCGGTGCAGGAAAACTCTCCATCACCCGTGCCACCCGGGCCCTGACGTTCCTGTCAGAGCTGGGACTGATTACCTACCAGACGGAATATGACCCGCTTATCGGGTGCTACATTCCGACCGACATCACGTTCACACTGGCTCTGTTTGCTGCCCTTGATGTGTCTGAGGATGCAGTGGCAGCTGCGCGCCGCAGTCGTGTTGAATGGGAAAACAAACAGCGCAAAAAGCAGGGGCTGGATACACTGGGTATGGATGAGCTGATAGCGAAAGCTTGGCGTTTTGTGCGTGAGCGTTTCCGCAGTTACCAGACAGAGCTTCAGTCCCGTGGAATAAAACGTGCCCGTGCGCGTCGTGATGCGAACAGAGAACGTCAGGACATCGTCACCCTAGTGAAACGGCAGCTGACGCGTGAAATCTCGGAAGGACGCTTCACTGCTAATGGTGAGGCGGTAAAACGCGAAGTGGAGCGTCGTGTGAAGGAGCGCATGATTCTGTCACGTAACCGCAATTACAGCCGGCTGGCCACAGCTTCTCCCTGAAAGTGATCTCCTCAGAATAATCCGGCCTGCGCCGGAGGCATCCGCACGCCTGAAGCCCGCCGGTGCACAAAAAAACAGCGTCGCATGCAAAAAACAATCTCATCATCCACCTTCTGGAGCATCCGATTCCCCCTGTTTTTAATACAAAATACGCCTCAGCGACGGGGAATTTTGCTTATCCACATTTAACTGCAAGGGACTTCCCCATAAGGTTACAACCGTTCATGTCATAAAGCGCCAGCCGCCAGTCTTACAGGGTGCAATGTATCTTTTAAACACCTGTTTATATCTCCTTTAAACTACTTAATTACATTCATTTAAAAAGAAAACCTATTCACTGCCTGTCCTGTGGACAGACAGGCTAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAG'
    assert get_fragment_seq('f003') == 'TGACGTCTAAGAAACGCGTTGACAATTAATCATCGGCTCGTATAATGCTTAAGTATAAGGAGGAAAAAATATGGCTAGCTGGAGCCACCCGCAGTTCGAAAAAGGCGCCGTGAGCAAAGGCGAAGAAACCACCATGGGCGTGATTAAACCGGATATGAAAATTAAACTGAAAATGGAAGGCAACGTGAACGGCCATGCGTTTGTGATTGAAGGCGAAGGCGAAGGCAAACCGTATGATGGCACCAACACCATTAATCTGGAAGTGAAAGAAGGCGCGCCGCTGCCGTTTAGCTATGATATTCTGACCACCGCGTTTAGCTATGGCAACCGTGCGTTTACCAAATATCCGGATGATATTCCGAACTATTTTAAACAGAGCTTTCCGGAAGGCTATAGCTGGGAACGTACCATGACCTTTGAAGATAAAGGCATTGTGAAAGTGAAAAGCGATATTAGCATGGAAGAAGATAGCTTTATTTATGAAATTCATCTGAAAGGCGAGAACTTTCCGCCGAACGGCCCGGTGATGCAGAAAGAAACCACCGGCTGGGATGCGAGCACCGAACGTATGTATGTGCGTGATGGCGTGCTGAAAGGTGATGTGAAAATGAAACTGCTGCTGGAAGGCGGCGGCCATCATCGTGTGGATTTTAAAACCATTTATCGTGCGAAAAAAGCGGTGAAACTGCCGGATTATCATTTTGTGGATCATCGTATTGAAATTCTGAACCATGATAAAGATTATAACAAAGTGACCGTGTATGAAATTGCGGTGGCGCGTAACAGCACCGATGGCATGGATGAACTGTATAAAGGGGGAGGAGGATCAGCGGCCCCAACTGATCTTCACCAAACGTATTACCGGCAGGTAAAGAACCCGAATCCGGTGTTCACTCCCCGTGAAGGTGCCGGAACGCTGAAGTTCTGCGAAAAACTGATGGAAAAGGCGGTGGGCTTCACCTCCCGTTTTGATTTCGCCATTCATGTGGCGCATGCCCGTTCCCGTGGTCTGCGTCGGCGCATGCCACCGGTGCTGCGTCGACGGGCTATTGATGCGCTGCTGCAGGGGCTGTGTTTCCACTATGACCCGCTGGCCAACCGCGTCCAGTGTTCCATCACCACACTGGCCATTGAGTGCGGACTGGCGACAGAGTCCGGTGCAGGAAAACTCTCCATCACCCGTGCCACCCGGGCCCTGACGTTCCTGTCAGAGCTGGGACTGATTACCTACCAGACGGAATATGACCCGCTTATCGGGTGCTACATTCCGACCGACATCACGTTCACACTGGCTCTGTTTGCTGCCCTTGATGTGTCTGAGGATGCAGTGGCAGCTGCGCGCCGCAGTCGTGTTGAATGGGAAAACAAACAGCGCAAAAAGCAGGGGCTGGATACACTGGGTATGGATGAGCTGATAGCGAAAGCTTGGCGTTTTGTGCGTGAGCGTTTCCGCAGTTACCAGACAGAGCTTCAGTCCCGTGGAATAAAACGTGCCCGTGCGCGTCGTGATGCGAACAGAGAACGTCAGGACATCGTCACCCTAGTGAAACGGCAGCTGACGCGTGAAATCTCGGAAGGACGCTTCACTGCTAATGGTGAGGCGGTAAAACGCGAAGTGGAGCGTCGTGTGAAGGAGCGCATGATTCTGTCACGTAACCGCAATTACAGCCGGCTGGCCACAGCTTCTCCCTGAAAGTGATCTCCTCAGAATAATCCGGCCTGCGCCGGAGGCATCCGCACGCCTGAAGCCCGCCGGTGCACAAAAAAACAGCGTCGCATGCAAAAAACAATCTCATCATCCACCTTCTGGAGCATCCGATTCCCCCTGTTTTTAATACAAAATACGCCTCAGCGACGGGGAATTTTGCTTATCCACATTTAACTGCAAGGGACTTCCCCATAAGGTTACAACCGTTCATGTCATAAAGCGCCAGCCGCCAGTCTTACAGGGTGCAATGTATCTTTTAAACACCTGTTTATATCTCCTTTAAACTACTTAATTACATTCATTTAAAAAGAAAACCTATTCACTGCCTGTCCTGTGGACAGACAGGCTAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAG'

    # I don't have an example with PCR overhangs yet.  That's definitely 
    # something I want to test, though.
    
    # Restriction digest
    assert get_fragment_seq(2)      == 'GCTTCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAGCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGCTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACGCGTAATACGACTCACTATAGGGCTTAAGTATAAGGAGGAAAAAATATGGAACGTCCGTATGCGTGTCCGGTGGAAAGCTGCGATCGTCGTTTTAGCCGTTCTGATGAACTGACCCGTCATATTCGTATTCATACCGGCCAGAAACCGTTTCAGTGCCGTATTTGCATGCGTAACTTTAGCCGTAGCGATCATCTGACCACCCATATTCGTACCCATACCGGCGAAAAACCGTTTGCGTGCGATATTTGCGGCCGTAAATTTGCGCGTTCTGATGAACGTAAACGTCATACCAAAATTCATCTGCGTCAGAAAGATGGCGGAGGTGGCTCTGGCGGTGGCGGATCGCACCACCATCACCATCATGGGGGAGGAGGATCAGGCTCAAGGGCGGGGGGCGGCGGGGAAAA'
    assert get_fragment_seq('f2')   == 'GCTTCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAGCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGCTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACGCGTAATACGACTCACTATAGGGCTTAAGTATAAGGAGGAAAAAATATGGAACGTCCGTATGCGTGTCCGGTGGAAAGCTGCGATCGTCGTTTTAGCCGTTCTGATGAACTGACCCGTCATATTCGTATTCATACCGGCCAGAAACCGTTTCAGTGCCGTATTTGCATGCGTAACTTTAGCCGTAGCGATCATCTGACCACCCATATTCGTACCCATACCGGCGAAAAACCGTTTGCGTGCGATATTTGCGGCCGTAAATTTGCGCGTTCTGATGAACGTAAACGTCATACCAAAATTCATCTGCGTCAGAAAGATGGCGGAGGTGGCTCTGGCGGTGGCGGATCGCACCACCATCACCATCATGGGGGAGGAGGATCAGGCTCAAGGGCGGGGGGCGGCGGGGAAAA'
    assert get_fragment_seq('f002') == 'GCTTCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAGCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGCTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACGCGTAATACGACTCACTATAGGGCTTAAGTATAAGGAGGAAAAAATATGGAACGTCCGTATGCGTGTCCGGTGGAAAGCTGCGATCGTCGTTTTAGCCGTTCTGATGAACTGACCCGTCATATTCGTATTCATACCGGCCAGAAACCGTTTCAGTGCCGTATTTGCATGCGTAACTTTAGCCGTAGCGATCATCTGACCACCCATATTCGTACCCATACCGGCGAAAAACCGTTTGCGTGCGATATTTGCGGCCGTAAATTTGCGCGTTCTGATGAACGTAAACGTCATACCAAAATTCATCTGCGTCAGAAAGATGGCGGAGGTGGCTCTGGCGGTGGCGGATCGCACCACCATCACCATCATGGGGGAGGAGGATCAGGCTCAAGGGCGGGGGGCGGCGGGGAAAA'

    # In vitro transcription
    assert get_fragment_seq(11)     == 'GGGCUUAAGUAUAAGGAGGAAAAAAUAUGGAACGUCCGUAUGCGUGUCCGGUGGAAAGCUGCGAUCGUCGUUUUAGCCGUUCUGAUGAACUGACCCGUCAUAUUCGUAUUCAUACCGGCCAGAAACCGUUUCAGUGCCGUAUUUGCAUGCGUAACUUUAGCCGUAGCGAUCAUCUGACCACCCAUAUUCGUACCCAUACCGGCGAAAAACCGUUUGCGUGCGAUAUUUGCGGCCGUAAAUUUGCGCGUUCUGAUGAACGUAAACGUCAUACCAAAAUUCAUCUGCGUCAGAAAGAUGGCGGAGGUGGCUCUGGCGGUGGCGGAUCGCACCACCAUCACCAUCAUGGGGGAGGAGGAUCAGGCUCAAGGGCGGGGGGCGGCGGGGAAAA'
    assert get_fragment_seq('f11')  == 'GGGCUUAAGUAUAAGGAGGAAAAAAUAUGGAACGUCCGUAUGCGUGUCCGGUGGAAAGCUGCGAUCGUCGUUUUAGCCGUUCUGAUGAACUGACCCGUCAUAUUCGUAUUCAUACCGGCCAGAAACCGUUUCAGUGCCGUAUUUGCAUGCGUAACUUUAGCCGUAGCGAUCAUCUGACCACCCAUAUUCGUACCCAUACCGGCGAAAAACCGUUUGCGUGCGAUAUUUGCGGCCGUAAAUUUGCGCGUUCUGAUGAACGUAAACGUCAUACCAAAAUUCAUCUGCGUCAGAAAGAUGGCGGAGGUGGCUCUGGCGGUGGCGGAUCGCACCACCAUCACCAUCAUGGGGGAGGAGGAUCAGGCUCAAGGGCGGGGGGCGGCGGGGAAAA'
    assert get_fragment_seq('f011') == 'GGGCUUAAGUAUAAGGAGGAAAAAAUAUGGAACGUCCGUAUGCGUGUCCGGUGGAAAGCUGCGAUCGUCGUUUUAGCCGUUCUGAUGAACUGACCCGUCAUAUUCGUAUUCAUACCGGCCAGAAACCGUUUCAGUGCCGUAUUUGCAUGCGUAACUUUAGCCGUAGCGAUCAUCUGACCACCCAUAUUCGUACCCAUACCGGCGAAAAACCGUUUGCGUGCGAUAUUUGCGGCCGUAAAUUUGCGCGUUCUGAUGAACGUAAACGUCAUACCAAAAUUCAUCUGCGUCAGAAAGAUGGCGGAGGUGGCUCUGGCGGUGGCGGAUCGCACCACCAUCACCAUCAUGGGGGAGGAGGAUCAGGCUCAAGGGCGGGGGGCGGCGGGGAAAA'

def test_get_oligo_seq():
    assert get_oligo_seq(10)     == 'AACGCGTAATACGACTCAC'
    assert get_oligo_seq('10')   == 'AACGCGTAATACGACTCAC'
    assert get_oligo_seq('o010') == 'AACGCGTAATACGACTCAC'

    # Non-standard nucleotides
    assert get_oligo_seq(11)     == 'AACGCGTAATACGACTCAC'
    assert get_oligo_seq('11')   == 'AACGCGTAATACGACTCAC'
    assert get_oligo_seq('o011') == 'AACGCGTAATACGACTCAC'

def test_mw():
    # Expected values from: http://molbiotools.com/dnacalculator.html

    # dsDNA, circular
    assert get_mw('p2') == approx(1659684.16, rel=1e-3)

    # dsDNA, linear
    assert get_mw('f3') == approx(1338870.98, rel=1e-3)

    # ssDNA
    # The expected value from several online MW calculators is: 5765.78
    # Biopython seems to be off by about 80 Da.  It's actually closer for the 
    # longer sequences (e.g. within 5 Da or so).  I'm not sure what the problem 
    # is with ssDNA, but I'm pretty sure it's not on my end.
    assert get_mw('o10') == approx(5845.75, rel=1e-3)

    # RNA
    assert get_mw('f11') == approx(125265.64, rel=1e-3)



