#!/usr/bin/env python3
import pyrosetta

from pyrosetta import Vector1
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
        #'-out:levels core.pack.rotamer_set.RotamerSet_.extra_rotamers:400',
]

pyrosetta.init(' '.join(flags))

pose = pose_from_pdb('1aay.pdb')
sfxn = core.scoring.ScoreFunctionFactory.create_score_function('ref2015')

protein_sele = ResiduePropertySelector(core.chemical.PROTEIN)
dna_sele = ResiduePropertySelector(core.chemical.DNA)
interface_sele = InterGroupInterfaceByVectorSelector(protein_sele, dna_sele)
protein_interface_sele = AndResidueSelector(interface_sele, protein_sele)

ex = core.pack.task.operation.ExtraRotamersGenericRLT()
ex.ex1(True)
ex.ex2(True)
ex.extrachi_cutoff(1)

tf = core.pack.task.TaskFactory()
tf.push_back(core.pack.task.operation.RestrictToRepacking())
tf.push_back(core.pack.task.operation.OperateOnResidueSubset(
    ex, protein_interface_sele))

task = tf.create_task_and_apply_taskoperations(pose)
rotsets = core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(pose)
ig = core.pack.interaction_graph.AnnealableGraphBase()

basic.options.set_string_vector_option('out:levels', Vector1(['core.pack.rotamer_set.RotamerSet_.extra_rotamers:500']))

core.pack.pack_rotamers_setup(pose, sfxn, task, rotsets, ig)

print('*' * 80)
print('*' * 80)
print('*' * 80)

basic.options.set_string_vector_option('out:levels', Vector1(['core.pack.rotamer_set.RotamerSet_.extra_rotamers:100']))
#basic.options.set_string_vector_option('out:mute', Vector1(['core.pack.rotamer_set.RotamerSet_.extra_rotamers']))

core.pack.pack_rotamers_setup(pose, sfxn, task, rotsets, ig)

