#!/usr/bin/env python3

"""\
Relax the given nucleoprotein in the Rosetta score function, while keeping any 
backbone movements within a range consistent with the structure's B-factors.

Usage:
    dnp_relax_b <xtal> [options]

Options:
    -p --percent-within PERCENT  [default: 90]
        The percent of atoms that should have remain within the average mean 
        squared displacement defined by their B-factors.  This option provides 
        a knob (hopefully one that is meaningful and intuitive) to control how 
        much the model moves during relaxation.  Higher values will produce 
        less movement, as fewer atoms are allowed to move significantly, while 
        lower value will produce more movement.  

    -t --convergence-tolerance PERCENT  [default: 5]
        How closely the objective (see --percent-within) needs to be satisfied 
        before the optimization converges.  The defaults specify that 85-95% of 
        the atoms in the model must be consistent with their B-factors.

    -r --restraint-range LOW,HIGH  [default: 0.5,5.0]
        The tightest and loosest restraints that will be considered, 
        respectively.  The values refer to the width of a harmonic well that 
        will be applied to each backbone atom, and should be separated by a 
        comma.  You can help the optimization converge faster by specifying a 
        narrower range, but if the actual tightness needed to achieve the given 
        percentage of consistent atoms is outside this range, the simulation 
        will fail and you will need to manually specify a larger range.

    -o --output DIRECTORY  [default: relax_?]
        The directory to put the output files in.  Any question marks (?) in 
        the given name will be replaced with the name of the given crystal 
        structure, without the extension.

    -f --force
        If the specified output directory already exists, overwrite it.  By 
        default, nothing will be overwritten and the program will exit with a 
        message explaining the error.

Before using Rosetta to perform protein design, the input structure must be 
relaxed in the Rosetta score function.  The goal of the relaxation step is to 
produce a model that scores as low as possible while remaining as similar as 
possible to the input structure.  These two goals are contrary to each other: 
the more a model is allowed to move, the lower scores it can achieve.  Thus it 
is necessary to define how much movement is acceptable.

Fortunately, this exact information is encoded in crystallographic B-factors.  
Each B-factor is related to the average coordinate error for a particular atom 
by the following equation (where <u²> if the average squared coordinate error):

    B = sqrt(<u²> / 8π²)

This script places harmonic restraints on the coordinates of each backbone 
atom, then optimizes the tightness of those restraints to keep the amount of 
movement consistent with the B-factors of the input structure.  More 
specifically, after the structure is relaxed, the squared displacement between 
the initial and relaxed coordinates of each atom is calculated.  This 
measurement is directly compared to the mean square displacement calculated 
from the B-factor for the same atom (if that atom has a B-factor).  The 
tightness of the restraints (i.e. the width of the harmonic well) is then 
optimized such that the 90% (by default, see --percent-within) of the backbone 
and base-pair atoms moved less than the displacement encoded in the B-factor.

Note that base-pair atoms are included in the B-factor comparisons, even though 
they are not restrained.  This is because we don't want the base-pairs to move 
significantly during relaxation, even though they are allowed to pack.  In 
contrast, amino-acid sidechains are also allowed to pack, but are expected to 
possibly adopt new rotamers, and are therefore not included in the B-factor 
comparisons.

Note also that optimizing restraints in this manner requires that the structure 
be relaxed ~10 times.  As each individual relaxation is expensive, be prepared 
for this protocol to take a long time.  You can reduce the number of times the 
structure must be relaxed by specifying a looser convergence threshold 
(--convergence-threshold) or a narrower restraint range (--restraint-range).
"""

import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *

import numpy as np
import pandas as pd
from math import pi

import sys, logging, shutil, json, stat
from pathlib import Path
from appdirs import AppDirs

logger = logging.getLogger('dnp_relax_b')

def main():
    import docopt

    try:
        args = docopt.docopt(__doc__)
        work = Workspace(
                args['<xtal>'],
                output_dir=args['--output'],
                overwrite=args['--force'],
        )

        init_rosetta()
        optimize_restraints(
                work,
                target_percentile=float(args['--percent-within']),
                interval=tuple(map(float, args['--restraint-range'].split(','))),
                tolerance=float(args['--convergence-tolerance']) / 100,
        )

    except WorkspaceExists as err:
        logger.error(err)
    except KeyboardInterrupt:
        print()

def optimize_restraints(work, *,
        target_percentile=90, interval=(0.5, 5.0), tolerance=0.01):
    # This is a 1D scalar root-finding problem.  The two primary algorithms 
    # provided by scipy for this kind of optimization are:
    #
    # Secant method:
    #   - Pros: Isn't bracketed and doesn't require derivative.
    #   - Cons: Can bounce around, and not robust if the function is "bouncier" 
    #     than expected.  Also, the algorithm will sample negative numbers, so 
    #     the objective function needs to convert those numbers into restraint 
    #     parameters that remain >0.
    #
    # Brentq method:
    #   - Pros: Best performance (i.e. fewest function evaluations).
    #   - Cons: Requires bracket, and fails if the zero falls outside that 
    #     bracket.
    #   
    # The Brentq algorithm seems to be more robust overall, and the range can 
    # be set via an argument.

    from scipy.optimize import brentq

    iteration = 0
    initial_pose = pose_from_pdb(str(work.initial_path))

    def objective(cst_stdev):
        nonlocal iteration
        iteration += 1

        # The tightness of the coordinate restraints can only be set via the 
        # command-line.  Welcome to Rosetta...
        basic.options.set_real_option('relax:coord_cst_stdev', cst_stdev)

        # Relax the pose with the given restraint.
        relaxed_pose = core.pose.Pose(initial_pose)
        sfxn = relax_pose(relaxed_pose)

        # Determine how well the relaxed structure agrees with the B-factors.
        df = calc_atoms_within_b(
                initial_pose, relaxed_pose,
                target_percentile=target_percentile,
        )
        objective = np.percentile(df.dist2_diff, target_percentile)

        # Record/log this iteration.
        initial_score = sfxn(initial_pose)
        relaxed_score = sfxn(relaxed_pose)
        score_diff = relaxed_score - initial_score

        work.record(
                iteration, 
                cst_stdev,
                relaxed_pose,
                objective=objective,
                score_reu=relaxed_score,
                score_diff_reu=score_diff,
        )

        logger.info(f"Iteration:                        {iteration}")
        logger.info(f"Tried:                            {cst_stdev}")
        logger.info(f"Mean squared distance (actual):   {df.dist2_actual.mean()}")
        logger.info(f"Mean squared distance (B-factor): {df.dist2_b_factor.mean()}")
        logger.info(f"Objective:                        {objective}")
        logger.info(f"Score (REU):                      {relaxed_score}")
        logger.info(f"Score improvement (REU):          {score_diff}")

        return objective

    x, results = brentq(
            objective,
            interval[0], interval[1],
            xtol=tolerance,
            full_output=True,
    )

    work.finalize(x)
    
    logger.info(f"Converged:      {results.converged}")
    logger.info(f"Root:           {results.root}")
    logger.info(f"Iterations:     {results.iterations}")
    logger.info(f"Function calls: {results.function_calls}")

def calc_atoms_within_b(initial_pose, relaxed_pose, target_percentile=90):
    rows = []

    for i in range(1, initial_pose.size() + 1):
        residue = initial_pose.residue(i)

        for j in range(1, residue.natoms() + 1):
            id = core.id.AtomID(j, i)
            v1 = initial_pose.xyz(id)
            v2 = relaxed_pose.xyz(id)
            d = v1.distance(v2)

            rows.append(dict(
                residue=i,
                atom=j,
                distance = v1.distance(v2),
                b_factor = initial_pose.pdb_info().bfactor(i, j),
                is_backbone = j <= residue.last_backbone_atom(),
                is_dna = residue.is_DNA(),
            ))

    df = pd.DataFrame(rows)
    df = df[df.b_factor > 0]
    df = df[df.is_backbone | df.is_dna]
    df['dist2_actual'] = df.distance**2
    df['dist2_b_factor'] = df.b_factor / (8*pi**2)
    df['dist2_diff'] = df.dist2_actual - df.dist2_b_factor
    return df


def init_rosetta():
    flags = [
            '-constant_seed',
            '-dna:specificity:exclude_dna_dna off',
            '-relax:dna_move on',  
            '-relax:coord_cst_stdev 0.5',  # default: 0.5
            '-relax:constrain_relax_to_start_coords on',
            '-relax:ramp_constraints off',
    ]
    pyrosetta.init(' '.join(flags), set_logging_handler='logging')

def relax_pose(pose):
    sfxn = load_score_function()
    tf = load_task_factory()
    mm = load_move_map()
    relax = load_fast_relax(sfxn, tf, mm)

    relax.apply(pose)
    return sfxn

def load_score_function():
    return core.scoring.ScoreFunctionFactory.create_score_function('ref2015')

def load_extra_rotamers():
    ex = core.pack.task.operation.ExtraRotamersGeneric()

    # Use extra rotamers for all χ angles: Arg and Lys are important for DNA 
    # interfaces, so we want to sample those sidechains finely.
    #ex.ex1(True)
    #ex.ex2(True)
    #ex.ex3(True)
    #ex.ex4(True)

    # DNA rotamers
    ex.exdna_sample_level(core.pack.task.NO_EXTRA_CHI_SAMPLES) # 1 rot/nt

    return ex
    
def load_task_factory():
    tf = core.pack.task.TaskFactory()
    # FastRelax automatically adds the `IncludeCurrent` task operation, so 
    # adding it here is redundant, but hopefully helpful in terms of clarity.
    tf.push_back(core.pack.task.operation.IncludeCurrent())
    tf.push_back(core.pack.task.operation.RestrictToRepacking())
    tf.push_back(load_extra_rotamers())

    return tf

def load_move_map():
    mm = core.kinematics.MoveMap()
    mm.set_bb(True)
    mm.set_chi(True)

    return mm

def load_fast_relax(sfxn, tf, mm):
    relax = protocols.relax.FastRelax(sfxn)
    relax.set_task_factory(tf)
    relax.set_movemap(mm)

    # Both of the following must be set, or no restraints will be applied.
    relax.constrain_coords(True)
    relax.constrain_relax_to_start_coords(True)

    return relax


def get_expected_error(pose):
    atoms = atoms_from_pose(pose)

    # Discard atoms without B-factors.
    atoms = atoms[atoms.b_factor > 0]

    return atoms.rmsd.median()

def atoms_from_pose(pose):
    info = pose.pdb_info()
    records = []

    for i in range(1, info.nres()+1):
        for j in range(1, info.natoms(i)+1):
            records.append(dict(
                residue=i,
                atom=j,
                chain=info.chain(i),
                b_factor=info.bfactor(i, j),
            ))

    atoms = pd.DataFrame(records)

    # Calculate expected coordinate error from the B-factors.
    atoms['rmsd'] = rmsds_from_b_factors(atoms.b_factor)

    return atoms

def rmsds_from_b_factors(bs):
    return np.sqrt(bs / (8*pi**2))
    
class Workspace:

    def __init__(self, xtal_path, output_dir, overwrite=False):
        xtal_path = Path(xtal_path)
        self.root = Path(output_dir.replace('?', xtal_path.stem)).resolve()

        if self.root.exists():
            if overwrite:
                shutil.rmtree(self.root)
            else:
                raise WorkspaceExists(self.root)

        self.root.mkdir()

        # Copy the command-line into the workspace.
        run_path = self.root / 'run.sh'

        with run_path.open('w') as f:
            f.write(f'''\
#!/usr/bin/env sh
{' '.join(sys.argv)}
''')

        run_mod = run_path.stat().st_mode
        run_path.chmod(run_mod | stat.S_IEXEC)

        # Copy the input structure into the workspace.  I decided to copy 
        # rather than link because copies are more robust and therefore a 
        # better record of what the simulation actually did.  Furthermore, the 
        # symlink shouldn't be necessary to figure out which input was used, 
        # because the workspace should be given a descriptive name.
        self.initial_path = self.root / ('initial' + ''.join(xtal_path.suffixes))
        shutil.copy(xtal_path, self.initial_path)

        # Track every model generated, so we can link to the best one.
        self.models = {}

        # Setup a JSON file that will describe the models generated during 
        # optimization.
        self.meta = self.root / 'optimization.json'
        with self.meta.open('w') as f:
            json.dump({}, f)

        # Log everything to the workspace.
        log = logging.getLogger()
        log.setLevel('DEBUG')

        formatter = logging.Formatter('{asctime}\t{name}\t{levelname}\t{message}', style='{')
        handlers = [
                logging.StreamHandler(),
                logging.FileHandler(self.root / 'log'),
        ]

        for handler in log.handlers[:]:
            log.removeHandler(handler)

        for handler in handlers:
            handler.setFormatter(formatter)
            log.addHandler(handler)

    def record(self, iteration, x, pose, **kwargs):
        # Write the PDB.
        path = self.root / f'iteration_{iteration}.pdb'
        pose.dump_pdb(str(path))

        # Keep track of the model generated by this restraint.
        self.models[x] = path

        # Add to the JSON log.
        with self.meta.open() as f:
            meta = json.load(f)

        meta[iteration] = {
                'path': str(path),
                'cst_stdev': x,
                **kwargs,
        }

        with self.meta.open('w') as f:
            json.dump(meta, f)

    def finalize(self, x):
        target = self.models[x]
        symlink = self.root / ('relaxed' + ''.join(target.suffixes))
        symlink.symlink_to(target.relative_to(self.root))


class WorkspaceExists(Exception):

    def __init__(self, path):
        super().__init__(f"'{path}' already exists.  Use `-f` to overwrite it.")
