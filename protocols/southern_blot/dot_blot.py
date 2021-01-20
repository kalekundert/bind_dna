#!/usr/bin/env python3
# vim: tw=49

"""\
Dot blot samples of DNA.

Usage:
    dot_blot <samples> [-n <n>] [-d <µL>] [-D <conc>] [-t <time>] [--saline]

Options:
    -n --num-spots <n>
        The number of DNA samples to spot.  By default, this is inferred from 
        the <samples> argument.

    -d --dna-volume <µL>        [default: 3]
        The volume of each sample of DNA to prepare.  The default is 3 µL, 
        enough for 2 µL to be spotted, but more may be needed if the DNA is 
        dilute or if multiple replicates will be spotted.

    -D --dna-conc <conc>
        The concentration of DNA to use.  The doesn't affect the amount of DNA
        prepared, it's just explicit about how much DNA should be used.

    -t --uv-time <time>         [default: 30−300 sec]
        How long to irradiate the membrane with UV light, to cross-link the
        DNA.  Only applicable for `--saline`.

    --saline
        Use the saline (rather than alkaline) transfer buffer.  The alkaline 
        buffer simplifies the protocol by removing the need for the UV cross-
        linking step, but the saline buffer may be useful for better mimicking 
        the conditions of a Southern blot.
"""

import docopt
import autoprop
from inform import plural
from stepwise import Protocol, MasterMix

@autoprop
class DotBlot:
    buffers = {
            'saline': """\
                Reagent   Stock  Volume  MM?
                =======  ======  ======  ===
                DNA              2.1 µL
                SSC         20x  0.9 µL  yes
            """,
            'alkaline': """\
                Reagent   Stock   Volume  MM?
                =======  ======  =======  ===
                DNA              1.65 µL
                NaOH        1 M  1.20 µL  yes
                EDTA     200 mM  0.15 µL  yes
            """,
    }

    def __init__(self):
        self.samples = None
        self.num_spots = None
        self.dna_volume_uL = None
        self.dna_conc = None
        self.uv_time = None
        self.transfer_type = 'alkaline'

    def get_reaction(self):
        rxn = MasterMix(self.buffers[self.transfer_type])
        rxn.extra_percent = 20
        rxn.extra_min_volume = '1 µL'

        if isinstance(self.samples, int):
            rxn.num_reactions = self.samples
        else:
            rxn['DNA'].name = ','.join(self.samples)
            rxn.num_reactions = len(self.samples)

        if x := self.num_spots:
            rxn.num_reactions = x
        if x := self.dna_volume_uL:
            rxn.hold_ratios.volume = x, 'µL'
        if x := self.dna_conc:
            rxn['DNA'].stock_conc = x

        if self.transfer_type == 'saline':
            rxn.show_master_mix = False

        return rxn

    def get_protocol(self):
        rxn = self.reaction
        p = Protocol()
        p += """\
Prepare a HyBond N+ membrane [1]:

- Draw a 50x50 mm grid using a dull pencil.
- Float in distilled water and allow to sink.
- Leave soaking for 10 min.
"""
        p += f"""\
Prepare {plural(rxn.num_reactions):# DNA sample/s} [1]:

{rxn}

- Incubate at 95°C for 10 min.
"""
        p += """\
Spot DNA [1]:

- Suspend the wetted membrane (e.g. over the top 
  of an open plastic box).
- Spot 2 µL of prepared DNA into each grid cell.
- Allow to dry.
"""
        if self.transfer_type == 'alkaline':
            p.steps[-1] += """\
- Rinse briefly in 2x SSC.
- Allow to dry.
"""
        if self.transfer_type == 'saline':
            p += f"""\
Denature and fix the DNA [1]:

- Place the membrane on sheets of Whatman 3MM 
  paper soaked in the following solutions:

  - 1.5 M NaCl, 0.5 M NaOH for 10 min.
  - 1.0 M NaCl, 0.5 M Tris HCl, pH=7.0 for 5 min.
  - Not soaked, until dry.

- Irradiate with 254 nm UV light for {self.uv_time}.
"""
        p.footnotes[1] = """\
Brown (1993). doi:10.1002/0471142727.mb0209bs21

The NaOH is to chemically hybridize the DNA to the
membrane.
"""
        return p

def int_or_none(x):
    return int(x) if x is not None else None

def float_or_none(x):
    return float(x) if x is not None else None

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    dot = DotBlot()

    try:
        dot.samples = int(args['<samples>'])
    except ValueError:
        dot.samples = args['<samples>'].split(',')

    dot.num_spots = int_or_none(args['--num-spots'])
    dot.dna_volume_uL = float_or_none(args['--dna-volume'])
    dot.dna_conc = args['--dna-conc']
    dot.uv_time = args['--uv-time']
    dot.transfer_type = 'saline' if args['--saline'] else 'alkaline'

    print(dot.protocol)
