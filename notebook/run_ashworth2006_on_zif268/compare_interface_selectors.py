#!/usr/bin/env python3

import pyrosetta
import pandas as pd
import os, sh

from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.residue_selector import (
        ResiduePropertySelector,
        InterGroupInterfaceByVectorSelector,
        OrResidueSelector,
        AndResidueSelector,
)
from pyrosetta.rosetta.protocols.dna import (
        DnaInterfaceSelector
)

def select_by_cb_vector(pose):
    # The Cα-Cβ vector will be <0,0,0> for everything, because the DNA doesn't 
    # have Cα.  See `core::select::util::cbeta_vector()`.  This means that 
    # nothing will be selected on the basis of angle.  A residue would be 
    # selected based on angle if the dot product between Cα1-Cβ1 and Cβ1-Cβ2 is 
    # greater than a threshold (0.25 by default), but the product will always 
    # be 0 if one of the operands is a null vector.

    # To be included based on distance, the rules are different for protein and 
    # DNA residues.  For protein, a sidechain atom needs to be within 
    # `nearby_atom_cut` of any atom in the other group.  For DNA, C1' (the 
    # sugar carbon supporting the base, defined as the "neighbor" atom, see 
    # database/chemical/...) needs to be within `nearby_atom_cut` of any atom 
    # in the other group.  So DNA is a little more restrictive, because it only 
    # counts distances from one atom, not all atoms.

    protein_sele = ResiduePropertySelector(core.chemical.PROTEIN)
    dna_sele = ResiduePropertySelector(core.chemical.DNA)
    interface_sele = InterGroupInterfaceByVectorSelector(protein_sele, dna_sele)
    protein_interface_sele = AndResidueSelector(interface_sele, protein_sele)

    return protein_interface_sele.apply(pose)
def select_by_arg_rotamers(pose):
    sele = DnaInterfaceSelector()
    return sele.apply(pose)

def run_pymol(pdb, **bools):
    args = [
            '-d', 'show sticks',
    ]
    colors = 'orange', 'yellow', 'white', 'cyan'

    if not os.fork():
        raise SystemExit

    for color, (name, bools_) in zip(colors, bools.items()):
        args += [
                '-d', pymol_sele_from_bools(name, bools_),
                '-d', f'color {color}, {name} and elem C',
        ]

    args += [
            '-d', 'select none',
    ]
    sh.pymol(pdb, *args)

def pymol_sele_from_bools(name, bools):
    resis = set()
    for i, x in enumerate(bools, 1):
        if x:
            resis.add(pose.pdb_info().number(i))

    sele = 'resi ' + '+'.join(str(x) for x in resis) if resis else 'none'
    return f'select {name}, {sele}'


if __name__ == '__main__':
    pyrosetta.init()

    pdb = '1aay.pdb'
    pose = pose_from_pdb(pdb)

    run_pymol(pdb,
            igibv=select_by_cb_vector(pose),
            argrot=select_by_arg_rotamers(pose),
    )



