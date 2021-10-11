#!/usr/bin/env python3

"""
Transfer protein to a PVDF membrane using an iBlot device

Usage:
    transfer_iblot [-p <program>] [-t <min>]

Options:
    -p --program <name>     [default: P3]
        Which program to run.  Each program specifies a voltage and a time.

    -t --time <min>
        How long to run the transfer.  If not specified, the default time for 
        the chosen protocol will be used.
"""

import stepwise
import docopt
from stepwise import pl, ul

programs = {
        'P0': ('20V,23V,25V', 7),
        'P1': ('25V', 6),
        'P2': ('23V', 6),
        'P3': ('20V', 7),
        'P4': ('15V', 7),
        'P5': ('10V', 7),
        'P6': ('7.5V', 3),
        'P7': ('5V', 3),
        'P8': ('20V,23V,25V', 7),
        'P9': ('20V,5V', 8),
}

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    program = args['--program']
    voltage, time_min = programs[program]
    if t := args['--time']:
        time_min = t

    if program not in programs:
        raise ValueError(f"unknown program: {program}")

    p = stepwise.Protocol()

    f1 = pl("""\
            https://tinyurl.com/yb5cj5ac
            https://tinyurl.com/t89e2n43
            """, """\
            Do not touch the membrane or gel with bare or 
            gloved hands. This may contaminate the gel or 
            membrane and interfere with further analysis. 
            If you need to adjust the membrane, always use 
            tweezers. 
            """, """\
            Sometimes there is green discoloration around 
            the sides of the membrane after the transfer.  
            This is due to copper ions from the transfer 
            stacks being carried by liquids and deposited 
            on the membrane.  These deposits do not 
            interfere with downstream processes.  The 
            stained regions can be cut away, but membrane 
            washing typically results in their removal.  
            To minimize this effect, shake excess water 
            off the filter paper and buffer from the gel 
            before placing each on the stack.
    """)

    p += pl(
            f"Transfer proteins to a PVDF membrane via iBlot{p.add_footnotes(f1)}:",
            ul(
                pl(
                    "Place the anode stack (in its plastic tray) in the iBlot machine.",
                    ul(
                        "Remove any bubbles underneath the membrane."
                    ),
                    br='\n',
                ),
            ),
            ul(
                pl(
                    "Place the gel on the membrane.",
                    ul(
                        "Rinse the gel with water beforehand.",
                        "Remove any bubbles.",
                    ),
                    br='\n',
                ),
            ),
            ul(
                "Soak included filter paper in water",
                pl(
                    "Place filter paper on gel.",
                    ul(
                        "Remove any bubbles.",
                    ),
                    br='\n',
                ),
            ),
            ul(
                pl(
                    "Place the cathode stack on the filter paper (copper side up, agarose side down, discard tray).",
                    ul(
                        "Remove any bubbles.",
                    ),
                    br='\n',
                ),
            ),
            ul(
                "Attach disposable sponge to the lid of the iBlot machine.  Make sure the metal band is in the upper right, aligned with the electrode.",
            ),
            ul(
                f"Run at {voltage} ({program}) for {time_min} min.",
            ),
            ul(
                "Transfer the membrane to water immediately after the run finishes.  If the membrane dries out, re-wet it in methanol and rinse with water before blocking.",
            ),
    )
    p.print()
