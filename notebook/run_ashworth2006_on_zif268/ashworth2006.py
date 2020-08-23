#!/usr/bin/env python3.6

from pyrosetta import init
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *

init('-dna:specificity:exclude_dna_dna off -constant_seed')

# Let's load Zif268:
pose = pose_from_pdb('1aay.pdb')

# Task Ops
# --------
# InitFromCommandLine: Nope
# IncludeCurrent: not sure if this is a good idea...maybe not for design?
# RestrictDesignToProteinDNAInterface: figure out exactly what this does.

# This should really be a ResidueSelector.

tf = core.pack.task.TaskFactory()
tf.push_back(
        protocols.dna.RestrictDesignToProteinDNAInterface()
)

# Score function
# --------------
# Thyme 2014 suggests using any DNA-optimized score function here.
# I'll of course also want to try ref2015.
sfxn = core.scoring.get_score_function();
print(sfxn.energy_method_options().exclude_DNA_DNA())
print(sfxn.energy_method_options().hbond_options().exclude_DNA_DNA())
sfxn.energy_method_options().exclude_DNA_DNA(False)
sfxn.energy_method_options().hbond_options().exclude_DNA_DNA(False)
print(sfxn.energy_method_options().exclude_DNA_DNA())
print(sfxn.energy_method_options().hbond_options().exclude_DNA_DNA())


raise SystemExit

# Movers
# ------
# DnaInterfacePacker: Gotta figure out what this does.
#
# This class is just really badly written...

# This pretty much just packs.
#
# If minimization is enabled (which it probably is, given the description of 
# the method in the paper), only sidechains to repackable residues are 
# minimized.  It's not clear to me if the original protocol repacks the DNA or 
# not.

mover = protocols.dna.DnaInterfacePacker()

mover.apply(pose)

pose.dump_pdb('1aay_pack_min.pdb')
