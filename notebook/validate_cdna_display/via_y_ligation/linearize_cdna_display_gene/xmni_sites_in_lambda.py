#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Restriction import XmnI
from Bio.SeqUtils import molecular_weight

records = SeqIO.parse('lambda_dna.fa', 'fasta')
record = next(records)

lambda_sites = len(XmnI.search(record.seq))
lambda_mw = molecular_weight(record.seq)

# 1U XmnI: Digest 1 µg λ DNA in 1h at 37°C
# Assume the plasmid is about the size of pKBK049 (2 kb), MW=1.5e6

plasmid_ug = 1
plasmid_mw = 1.5e6

xmni_units = plasmid_ug * (1 / plasmid_mw) * (lambda_mw / lambda_sites)

print(f"λ genome: {lambda_mw:.1f} Da, {lambda_sites} XmnI sites")
print(f"Use ≥{xmni_units:.2f} U XmnI to digest {plasmid_ug} µg ~2 kb plasmid DNA in 1h at 37°C")




