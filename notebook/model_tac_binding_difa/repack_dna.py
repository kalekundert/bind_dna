#!/usr/bin/env python3

import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *

from pyrosetta.rosetta.core.pack.task.operation import (
        IncludeCurrent,
        RestrictToRepacking,
        OperateOnResidueSubset,
        PreventRepackingRLT,
        ExtraRotamersGeneric,
)
from pyrosetta.rosetta.core.select.residue_selector import (
        ResiduePropertySelector,
)

pyrosetta.init('''\
        -constant_seed
        -jran 0
        -dna:specificity:exclude_dna_dna off
''')

pose = pose_from_pdb('dna_only/unrelaxed.pdb')
sfxn = core.scoring.ScoreFunctionFactory.create_score_function('ref2015')

ex = ExtraRotamersGeneric()
#ex.exdna_sample_level(core.pack.task.EX_THREE_THIRD_STEP_STDDEVS) # 500 rot/nt
#ex.exdna_sample_level(core.pack.task.EX_TWO_FULL_STEP_STDDEVS) # 50 rot/nt
#ex.exdna_sample_level(core.pack.task.EX_ONE_HALF_STEP_STDDEV) # 25 rot/nt
#ex.exdna_sample_level(core.pack.task.EX_ONE_STDDEV) # 10 rot/nt
#ex.exdna_sample_level(core.pack.task.NO_EXTRA_CHI_SAMPLES) # 1 rot/nt
ex.exdna_sample_level(core.pack.task.EX_SIX_QUARTER_STEP_STDDEVS) # DB

tf = core.pack.task.TaskFactory()
tf.push_back(IncludeCurrent())
tf.push_back(RestrictToRepacking())
tf.push_back(ex)

pack = protocols.minimization_packing.PackRotamersMover()
pack.score_function(sfxn)
pack.task_factory(tf)

pack.apply(pose)

pose.dump_pdb('dna_only/repacked.pdb')
