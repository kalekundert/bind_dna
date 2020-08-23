#!/usr/bin/env python3

import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *

import numpy as np
import pandas as pd

pyrosetta.init()


p1 = pose_from_pdb('1bna.pdb')
p2 = pose_from_pdb('1bna_redo.pdb')

rows = []


for i in range(1, p1.size() + 1):

    r1 = p1.residue(i)
    n = r1.last_backbone_atom()
    print(i, r1.last_backbone_atom())

    for atom in r1.atoms():
        atom.show()



    for j in range(1, p1.residue(i).natoms() + 1):
        id = core.id.AtomID(j, i)
        v1 = p1.xyz(id)
        v2 = p2.xyz(id)
        d = v1.distance(v2)

        rows.append(dict(
            residue=i,
            atom=j,
            distance=d,
            b_factor=p1.pdb_info().bfactor(i, j),
        ))

df = pd.DataFrame(rows)
df = df[df.b_factor > 0]
df['b_disp'] = np.sqrt(df.b_factor / (8*np.pi**2))
df['delta_d'] = df.b_disp - df.distance

df = df.set_index(['residue', 'atom'])

print(df.describe())

rmsd = np.sqrt((df.distance**2).mean())
print(rmsd)
