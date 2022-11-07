#!/usr/bin/env python3

from Bio import Restriction
from Bio.Seq import Seq

primers = {
        'SR008': 'GTATGTCGGCTCTCGTATCG',
        'SR013': 'GCAGCGTTTTAGCCTACAAG',
        'SR015': 'GAAGTGGTTTGCCTAAACGC',
        'SR022': 'GGGTTGTCTCCTCTGATAGC',
        'SR027': 'CGGGAGGAAGTCTTTAGACC',
        'SR031': 'CTAATATCCCTGAGCGACGG',
        'SR034': 'GGAAACAATAACCATCGGCG',
        'SR038': 'CCACGAGATAAGAGGATGGC',
        'SR041': 'GTATAAGATCAGCCGGACCC',
        'SR042': 'CCTTTAACAGGACATGCAGC',
        'SR045': 'CGGATCGAACTTAGGTAGCC',
        'SR054': 'GCTAAATAGAGGGAAGCCCC',
        'SR056': 'GCATAAAGTTGACAGGCCAG',
        'SR059': 'CGATAGAACGACCAGGTAGC',
        'SR065': 'CCCGAGGGGAGAAATATACC',
        'SR069': 'GGAAAACTAAGACAAGGCGC',
        'SR071': 'CCGGTTGTACCTATCGAGTG',
        'SR075': 'GTTGCATCTAAGCCAAGTGC',
        'SR078': 'CTATAGAATCCGGGCTGGTC',
        'SR080': 'GGGCACCGATTAAGAAATGC',
        'SR086': 'CACTCGATAGGTACAACCGG',
        'SR091': 'GTACTCAGAGATTGCCGGAG',
        'SR093': 'GCACGCAAAAGGACATAACC',
        'SR096': 'CAGACCTACGGATCTTAGCG',
        'SR098': 'CCAGAGCTTAGGGGACATAC',
        'SR113': 'GACCATGCAAGGAGAGGTAC',
        'SR119': 'CCGGGAGGAAGATATAGCAC',
        'SR120': 'CCGTGCGACAAGATTTCAAG',
        'SR122': 'CCGAGGGAACCATGATACAG',
        'SR123': 'GAAAAGTCCCAATGAGTGCC',
        'SR133': 'GTTCAGAGGTACGAACCCTC',
        'SR139': 'GATACATAGACTTGGCCCCG',
        'SR146': 'CGAACGCAAAAGTCCTCAAG',
        'SR151': 'CTAGGGGATGGTCCAATACG',
        'SR159': 'CTGCTAGGGGCTACTTATCG',
        'SR163': 'CTAGGGAACCAGGCTTAACG',
        'test': 'GAATTC',
}
enz = [
        Restriction.AgeI,
        Restriction.ApoI,
        Restriction.BamHI,
        Restriction.BbsI,
        Restriction.BclI,
        Restriction.BmtI,
        Restriction.BsaI,
        Restriction.BsiWI,
        Restriction.BsrGI,
        Restriction.BstEII,
        Restriction.BstZ17I,
        Restriction.DraIII,
        Restriction.EagI,
        Restriction.EcoRI,
        Restriction.EcoRV,
        Restriction.HindIII,
        Restriction.KpnI,
        Restriction.MfeI,
        Restriction.MluI,
        Restriction.NcoI,
        Restriction.NheI,
        Restriction.NotI,
        Restriction.NruI,
        Restriction.NsiI,
        Restriction.PstI,
        Restriction.PvuI,
        Restriction.PvuII,
        Restriction.SacI,
        Restriction.SalI,
        Restriction.SbfI,
        Restriction.ScaI,
        Restriction.SpeI,
        Restriction.SphI,
        Restriction.SspI,
        Restriction.StyI,
]

for name, primer in primers.items():
    seq = Seq(primer)
    hits = []

    for e in enz:
        if e.search(seq):
            hits.append(e)

    if hits:
        print(name, hits)
