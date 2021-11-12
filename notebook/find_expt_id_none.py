#!/usr/bin/env python3

import exmemo
from pathlib import Path

work = exmemo.Workspace.from_path(__file__)

for expt in work.iter_experiments():
    print(expt.id)
