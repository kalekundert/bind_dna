#!/usr/bin/env python3

import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *

# The default algorithm for setting the fold tree somehow isn't compatible with 
# DNA.  This produces a whole bunch of error messages, but since the DNA 
# doesn't move, it shouldn't cause any real problems.  Possible solutions:
# 
# - Mute the relevant tracer (core.kinematics.AtomTree).  This is appropriate 
#   because the error isn't really a problem.
# 
# - Trust the fold tree: ``lm.trust_fold_tree()``.  This does get rid of the 
#   error message, but it would incur a performance penalty (because more 
#   torsions would have to be updated on each move), and would introduce the 
#   possibility of atoms throughout the structure moving slightly.
#
# - Manually make a fold tree and trust it.  This is probably the "proper" 
#   thing to do, but it'd be the most effort.
#
# I'm just going to mute the tracer for now.  That seems like the most 
# proportionate response for this simple script.

pyrosetta.init('''\
        -constant_seed
        -jran 0
        -loops.frag_sizes 9 3 1
        -loops.frag_files frags/aa5hxy_09_05.200_v1_3 frags/aa5hxy_03_05.200_v1_3 none
        -dna:specificity:exclude_dna_dna off
        -mute core.kinematics.AtomTree
''')

pose = pose_from_pdb('005_fix_difA.pdb')
resi = pose.pdb_info().pdb2pose

# Chain B loop: Not finding fragments for some reason...

loops = protocols.loops.Loops()
loops.add_loop(resi('A', 87), resi('A', 106), resi('A', 97), 0, 1)
#loops.add_loop(resi('B', 87), resi('B', 106), resi('B', 97), 0, 1)

lm = protocols.loop_modeler.LoopModeler()
lm.setup_kic_with_fragments_config()
#lm.disable_centroid_stage()
lm.disable_fullatom_stage()
lm.set_loops(loops)
lm.apply(pose)

pose.dump_pdb('006_connect_domains.pdb')
