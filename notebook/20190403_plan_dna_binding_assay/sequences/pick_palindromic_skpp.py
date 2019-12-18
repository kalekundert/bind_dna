#!/usr/bin/env python3

"""\
Usage:
    pick_palindromic_skpp.py [<n>]

Arguments:
    <n>
        The number of primers to display.  Primers are sorted based on three 
        criteria: how palindromic they were to begin with, how near to 50% GC 
        content they are, and how long of a 3' GC clamp they have.
"""

import pandas as pd, re, docopt
from textdistance import levenshtein

skpp_csv = '../skpp/patent_20140045728_skpp.csv'

def complement(sequence):
    complements = str.maketrans('ACTGactg', 'TGACtgac')
    return sequence.translate(complements)

def reverse_complement(sequence):
    return complement(sequence[::-1])

def gc_content(sequence):
    return sum(x in 'GC' for x in sequence.upper()) / len(sequence)

def gc_clamp(sequence):
    m = re.search('[GC]*$', sequence)
    return m.end() - m.start()


def load_primers():
    df = pd.read_csv(skpp_csv, header=None, names=['name', 'primer'])
    df = make_palindromes(df)
    return df

def rank_primers(df):
    df['levenshtein'] = df.apply(
            lambda x: levenshtein(x.primer, x.parent),
            axis=1,
    )
    df['gc_content'] = df.apply(
            lambda x: gc_content(x.primer),
            axis=1,
    )
    df['gc_clamp'] = df.apply(
            lambda x: gc_clamp(x.primer),
            axis=1,
    )
    df['gc_balance'] = abs(df.gc_content - 0.5)

    df.sort_values(
            ['gc_balance', 'levenshtein', 'gc_clamp'],
            ascending=[True, True, False],
            inplace=True,
    )
    df.drop('gc_balance', axis='columns', inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df

def make_palindromes(df):
    def make_palindromes_5(row):
        i = len(row.primer) // 2
        return pd.Series({
            'parent': row.primer,
            'name': row['name'] + '-5',
            'primer': row.primer[:i] + reverse_complement(row.primer)[i:],
        })

    def make_palindromes_3(row):
        i = len(row.primer) // 2
        return pd.Series({
            'parent': row.primer,
            'name': row['name'] + '-3',
            'primer': reverse_complement(row.primer)[:i] + row.primer[i:] 
        })

    df = pd.concat([
            df.apply(make_palindromes_5, axis=1),
            df.apply(make_palindromes_3, axis=1),
        ],
    )
    return df.sort_values('name').reset_index(drop=True)


if __name__ == '__main__':
    import docopt
    args = docopt.docopt(__doc__)

    df = load_primers()
    df = rank_primers(df)

    print(df.head(int(args['<n>'] or 10)))




