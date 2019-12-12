#!/usr/bin/env python

import pyrosetta
import pandas as pd

from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.residue_selector import (
        ResiduePropertySelector,
        InterGroupInterfaceByVectorSelector,
        OrResidueSelector,
        AndResidueSelector,
)

flags = [
        '-constant_seed',

        # Score interactions between DNA atoms.
        '-dna:specificity:exclude_dna_dna off',

        # Allow DNA to move during FastRelax.  This is undocumented, and can 
        # only be set via command-line option.  It also kinda suggests that 
        # relaxing structures with DNA might not work very well...
	'-relax:dna_move on',  

        # Can't be set programatically.  Would need to instantiate 
        # AtomCoordinateCstMover and replicate everything FastRelaxMover does 
        # with that object.
        '-relax:coord_cst_stdev 0.5',  # Default: 0.5

        #'-relax:constrain_relax_to_start_coords on'

        '-relax:ramp_constraints off',

        #'-out:levels core.pack.rotamer_set.RotamerSet_:400',
        '-out:levels core.pack.rotamer_set.RotamerSet_.extra_rotamers:400',
]

pyrosetta.init(' '.join(flags))

def print_sfxn(sfxn, pose):
    print("Exclude DNA/DNA:        ", sfxn.energy_method_options().exclude_DNA_DNA())
    print("Exclude DNA/DNA H-bonds:", sfxn.energy_method_options().hbond_options().exclude_DNA_DNA())
    print("Score:                  ", sfxn(pose), "REU")

    s = std.ostringstream()
    pose.energies().show_total_headers(s)
    headers = s.str().split()

    s = std.ostringstream()
    pose.energies().show_totals(s)
    totals = [float(x) for x in s.str().split()]

    df = pd.DataFrame({'term': headers, 'unweighted': totals})
    print(df)
    print()

pose = pose_from_pdb('1aay.pdb')

#x = protocols.dna.DnaInterfaceFinder()
#x.determine_protein_interface(pose)
#
#n = x.protein_neighbors()

def pymol_from_sele(name, pose):
    bools = globals()[name].apply(pose)
    resis = set()
    for i, x in enumerate(bools, 1):
        if x:
            resis.add(pose.pdb_info().number(i))


    sele = 'resi ' + '+'.join(str(x) for x in resis) if resis else 'none'
    print(f'select {name}, {sele}')

protein_sele = ResiduePropertySelector(core.chemical.PROTEIN)
dna_sele = ResiduePropertySelector(core.chemical.DNA)
interface_sele = InterGroupInterfaceByVectorSelector(protein_sele, dna_sele)

# The Cα-Cβ vector will be <0,0,0> for everything, because the DNA doesn't have 
# Cα.  See `core::select::util::cbeta_vector()`.  This means that nothing will 
# be selected on the basis of angle.  A residue would be selected based on 
# angle if the dot product between Cα1-Cβ1 and Cβ1-Cβ2 is greater than a 
# threshold (0.25 by default), but the product will always be 0 if one of the 
# operands is a null vector.

#interface_sele.vector_dist_cut(0.0)
#interface_sele.vector_angle_cut(0.0)
#interface_sele.nearby_atom_cut(100.0)
#interface_sele.cb_dist_cut(100.0)

# To be included based on distance, the rules are different for protein and DNA 
# residues.  For protein, a sidechain atom needs to be within `nearby_atom_cut`
# of any atom in the other group.  For DNA, C1' (the sugar carbon supporting 
# the base, defined as the "neighbor" atom, see database/chemical/...) needs to 
# be within `nearby_atom_cut` of any atom in the other group.  So DNA is a 
# little more restrictive, because it only counts distances from one atom, not 
# all atoms.

protein_interface_sele = AndResidueSelector(interface_sele, protein_sele)

print('interface_sele.vector_dist_cut =', interface_sele.vector_dist_cut())
print('interface_sele.vector_angle_cut =', interface_sele.vector_angle_cut())
print('interface_sele.nearby_atom_cut =', interface_sele.nearby_atom_cut())
print('interface_sele.cb_dist_cut =', interface_sele.cb_dist_cut())

pymol_from_sele('protein_sele', pose)
pymol_from_sele('dna_sele', pose)
pymol_from_sele('interface_sele', pose)
pymol_from_sele('protein_interface_sele', pose)

raise SystemExit

#n_close = 0
#n_contact = 0
#n_total = 0
#
#
#for k, v in n.items():
#    n_close += v.close()
#    n_contact += v.contact()
#    n_total += 1 
#    print(f"{pose.residue(k).name1()}{k:<4} close={int(v.close())} contact={int(v.contact())}")
#
#print()
#print(f"close: {100 * n_close / n_total}%")
#print(f"contact: {100 * n_contact / n_total}%")

## Scoring

sfxn = core.scoring.ScoreFunctionFactory.create_score_function('ref2015')
#sfxn = get_score_function()
#print_sfxn(sfxn, pose)

## Packing

ex = core.pack.task.operation.ExtraRotamersGeneric()
#ex.exdna_sample_level(core.pack.task.EX_THREE_THIRD_STEP_STDDEVS) # 500 rot/nt
#ex.exdna_sample_level(core.pack.task.EX_TWO_FULL_STEP_STDDEVS) # 50 rot/nt
#ex.exdna_sample_level(core.pack.task.EX_ONE_HALF_STEP_STDDEV) # 25 rot/nt
#ex.exdna_sample_level(core.pack.task.EX_ONE_STDDEV) # 10 rot/nt
ex.exdna_sample_level(core.pack.task.NO_EXTRA_CHI_SAMPLES) # 1 rot/nt
ex.ex4(True)
#ex.exdna_sample_level(core.pack.task.EX_SIX_QUARTER_STEP_STDDEVS) # DB

tf = core.pack.task.TaskFactory()
# FastRelax adds IncludeCurrent automatically, so this isn't necessary for tha 
# that application.
tf.push_back(core.pack.task.operation.IncludeCurrent())
tf.push_back(core.pack.task.operation.RestrictToRepacking())
tf.push_back(ex)

pack = protocols.minimization_packing.PackRotamersMover()
#pack = protocols.dna.DnaInterfacePacker()
pack.score_function(sfxn)
pack.task_factory(tf)

pack.apply(pose)
raise SystemExit

#print_sfxn(sfxn, pose)

## Minimization

mm = core.kinematics.MoveMap()
mm.set_bb(True)
mm.set_chi(True)

min = protocols.minimization_packing.MinMover()
min.score_function(sfxn)
min.movemap(mm)

#min.apply(pose)
#print_sfxn(sfxn, pose)

## Relaxing

# If you want constraints on for the entire relax run, set ramp_constraints to 
# false (below) along with the other constraint flags.  Constraints can be 
# provided using either the cst_fa_file or cst_file options.  If both options 
# are used, priority will be given to the cst_fa_file constraints. Default 
# weights for the constraints are 0. Built in options include backbone 
# coordinate constraints, sidechain coordinate constraints and sidechain 
# pairwise constraints. 

# -constraints:cst_fa_file  <filename>      
#
#    Add constraints from the fullatom constraint file(*)
#
# -constraints:cst_fa_weight <weight>
# 
#   Weight to be used for the constraints in the cst_fa_file.
#
# -constraints:cst_file  <filename>
#
#    Add constraints from the constraint file(*)
#
# -constraints:cst_weight   <weight>
#
#    Weight to be used for the constraints in the cst_file.
#
# -relax:constrain_relax_to_start_coords
#
#    Add coordinate constraints to backbone heavy atoms, based on the input 
#    structure.
#
# -relax:constrain_relax_to_native_coords
#
#    Add coordinate constraints to backbone heavy atoms, based on the structure 
#    passed to -in:file:native.
#
# -relax:coord_constrain_sidechains
#
#    Also add coordinate constraints to sidechain heavy atoms (requires one of 
#    previous two options)
#
# -relax:coord_cst_stdev    <stdev>
#
#    Set the strength of coordinate constraints (smaller=tighter)
#
# -relax:coord_cst_width    <width>
#
#    If set, use flat-bottomed constraints instead of harmonic constraints, 
#    with a bottom width of <width>
#
# -relax:sc_cst_maxdist   <dist>
#
#    Add pairwise atom constraints to sidechain atoms up to dist apart from one 
#    another.
#
# -relax:ramp_constraints   false
#
#    When explicitly set to false, do not ramp down constraints (does not 
#    affect ramping in custom scripts). When true, constraints are ramped down 
#    during each simulated annealing cycle.
#
# -relax:dualspace true
#
#    Use the Dualspace protocol for dihedral and cartesian minimization as 
#    described by Conway et al. Do 3 FastRelax cycles of internal coordinate 
#    relax followed by two cycles of Cartesian relax - cart_bonded energy term 
#    is required, pro_close energy term should be turned off, and use of 
#    -relax::minimize_bond_angles is recommended.  The -nonideal flag can be 
#    used for this.


#relax = protocols.relax.FastRelax(sfxn, 'MonomerRelax2019')
#relax = protocols.relax.FastRelax(sfxn, 'my_fastrelax.txt')
relax = protocols.relax.FastRelax(sfxn)
relax.set_task_factory(tf)
relax.set_movemap(mm)

# Both of the following must be set, or no restraints will be applied.
relax.constrain_coords(True)
relax.constrain_relax_to_start_coords(True)

relax.ramp_down_constraints(False)
relax.cartesian(True)

relax.apply(pose)
print_sfxn(sfxn, pose)


pose.dump_pdb('1bna_repack.pdb')

