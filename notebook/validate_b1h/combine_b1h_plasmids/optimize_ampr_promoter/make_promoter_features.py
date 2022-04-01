#!/usr/bin/env python3

"""\
Annotate predicted promoters in the given plasmid, using predictions 
made by the Salis Lab Promoter calculator.

Usage:
    make_promoter_features.py <plasmid.dna> <promoters.csv> [-o <path>]

Arguments:
    <plasmid.dna>
        A Snapgene plasmid file that you wish to annotate.

    <promoters.csv>
        A CSV file containing promoter predictions.  This file can be 
        downloaded from the "De Novo DNA" web server.  The following columns 
        must be present:

        - TSS  [Transcriptional Start Site Position (nt)]
        - Tx_rate  [Transcription Initiation Rate (au)]
        - sequence [promoter sequence]

Options:
    -o --output <path>
        Where to write the modified Snapgene plasmid file.

Details:
    Predicted promoters are added to the plasmid in order of strength, as long 
    as they do not overlap with any other predicted promoters.  The process 
    stops once promoters overlapping all previously annotated promoters in the 
    file have been found.
"""

import pandas as pd
import autosnapgene as snap

from matplotlib.cm import get_cmap
from matplotlib.colors import to_hex

if __name__ == '__main__':
    import docopt
    args = docopt.docopt(__doc__)
    dna = snap.parse(args['<plasmid.dna>'])
    df = pd.read_csv(args['<promoters.csv>'])
    cm = get_cmap('viridis')

    df = df.sort_values(
            'Tx_rate  [Transcription Initiation Rate (au)]',
            ascending=False,
    )

    promoters = [x for x in dna.features if x.type == 'promoter']
    unlabeled_promoters = promoters[:]

    # Better way to avoid duplicates: group by -10 and -35, then pick the best 
    # score in each group.
    already_seen = {
            '+': set(),
            '-': set(),
    }

    max_rate = df.iloc[0]['Tx_rate  [Transcription Initiation Rate (au)]']
    min_rate = df.iloc[-1]['Tx_rate  [Transcription Initiation Rate (au)]']
    n = 0

    for i, row in df.iterrows():
        seq = row['sequence [promoter sequence]']
        rate = row['Tx_rate  [Transcription Initiation Rate (au)]']
        tss = int(row['TSS  [Transcriptional Start Site Position (nt)]'])
        itr = row['ITR [initial transcribed region sequence]']
        strand = row['strand  [Strand Orientation]']
        direction = {'+': 'forward', '-': 'backward'}[strand]

        start = tss - len(itr)
        end = start + len(seq)
        indices = range(start, end)

        rel_rate = (rate - min_rate) / (max_rate - min_rate)

        if already_seen[strand].intersection(indices):
            continue

        already_seen[strand].update(indices)

        feat = snap.Feature.from_segment(
                name=f'P{tss}',
                type='promoter',
                range=(start, end),
                directionality=direction,
                color=to_hex(cm(rel_rate)),
                qualifiers={
                    'note': f'{rate:.2f} au',
                },
        )
        feats = dna.add_feature(feat)
        n += 1

        unlabeled_promoters = [
                f for f in unlabeled_promoters
                if not (
                    set(range(f.range[0], f.range[1] + 1)).intersection(indices) and
                    f.directionality == direction
                )
        ]
        if not unlabeled_promoters:
            print(f"Added {n} promoters.")
            print(f"Worst promoter: P{tss}")
            print(f"  Rate: {rate:.2f} au")
            break

    dna.write(args['--output'])
