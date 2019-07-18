#!/usr/bin/env python3

"""\
Relax the given nucleoprotein in the Rosetta score function, while keeping any 
backbone movements within a range consistent with the structure's B-factors.

Usage:
    dbp_relax_b init [<workspace>] <xtal> [-fd]
    dbp_relax_b optimize <workspace> [-fP]
    dbp_relax_b relax <workspace> [<id>] [-f]
    dbp_relax_b finish <workspace> [-f]

Example:
    The first step is to create a new workspace (i.e. a directory which will 
    contain all the input and output files for the simulations):

    $ dbp_relax_b init 1aay.pdb

    The second step is to determine how strong the backbone restraints need to 
    be to produce the desired amount of movement.  Expect this step to 
    take ~1 day to complete:
    
    $ dbp_relax_b optimize relax_1aay

    The third step is to repeat the relaxation enough times (with the optimized 
    restraints) to mitigate any stochastic effects.  Use your cluster to run 
    the below command as many times as necessary for the scores to converge 
    (typically 10-100).  Each simulation may take a few hours:

    $ dbp_relax_b relax relax_1aay

    The last step is to identify the lowest scoring model, for use in future 
    design simulations:

    $ dbp_relax_b finish relax_1aay

Subcommands:
    init
        Create a "workspace", which is a directory containing all the relevant 
        input and output files.  The <xtal> argument is the path the PDB file 
        to relax.  If no workspace name is specified, the workspace will be 
        named `relax_<xtal>` (without the `.pdb` suffix).  A file called 
        `conf.toml` specifying and documenting all the default settings will 
        automatically be added to the new workspace.  You may edit this file to 
        change these settings.

    optimize
        Run a series of simulations to determine how strong the backbone 
        coordinate restraints should be to recapitulate the about of movement 
        described by the B-factors.  Because the optimization requires the 
        simulations to run serially, this step can take a long time.  

        Note that this command actually proceeds in two steps.  In the first 
        step, abbreviated FastRelax simulations are run to quickly identify 
        reasonable restraint weights.  The second step, full-length FastRelax 
        simulations are run to accurately identify the optimal restraint 
        weight.  You can use the --skip-preoptimize flag to skip the first 
        step.

    relax
        Run a single FastRelax simulation using the optimal restraint weight 
        determined previously.  It is an error to run this command before the 
        `optimize` command.  Typically you would run 10-100 of these 
        simulations using a cluster.

        Each `relax` simulation needs a unique id, or the results will 
        overwrite each other.  You can specify an id via the optional <id> 
        argument, otherwise value of the $SLURM_ARRAY_TASK_ID environment 
        variable will be used.

    finish
        Identify the relaxed model with the best score.  A symlink to this 
        model named `relaxed.pdb` will be created in the workspace directory.

Options:
    -P --skip-preoptimize
        Skip the first, approximate optimization when performing the `optimize` 
        command.  This is typically necessary if the root found by the first 
        optimization turns out to be far from the real root, in which case the 
        second optimization will fail.  In this case, you will need to either 
        increase the `optimize.preoptimize-delta` setting in `conf.toml`, or 
        delete the `preoptimize` directory in the workspace and manually set
        the `optimize.restraint-range` setting in `conf.toml` such that the 
        real root is within the given range.

    -f --force
        If the specified workspace already exists, overwrite it.  By default, 
        nothing will be overwritten and the program will exit with a message 
        explaining the error.

    -d --debug
        Indicate that the workspace is meant to be used for debugging.  All 
        simulations will be run with drastically reduced quality in the 
        interest of completing quickly.  Never use the results of these 
        simulations in downstream steps.

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

The optimize command places harmonic restraints on the 
coordinates of each backbone atom, then optimizes the tightness of those 
restraints to keep the amount of movement consistent with the B-factors of the 
input structure.  More specifically, after the structure is relaxed, the 
squared displacement between the initial and relaxed coordinates of each atom 
is calculated.  This measurement is directly compared to the mean square 
displacement calculated from the B-factor for the same atom (if that atom has a 
B-factor).  The tightness of the restraints (i.e. the width of the harmonic 
well) is then optimized such that the 90% (by default, see the `percent-within` 
setting) of the backbone and base-pair atoms moved less than the displacement 
encoded in the B-factor.

Note that base-pair atoms are included in the B-factor comparisons, even though 
they are not restrained.  This is because we don't want the base-pairs to move 
significantly during relaxation, even though they are allowed to pack.  In 
contrast, amino-acid sidechains are also allowed to pack, but are expected to 
possibly adopt new rotamers, and are therefore not included in the B-factor 
comparisons.

Note also that optimizing restraints in this manner requires that the structure 
be relaxed ~10 times.  As each individual relaxation is expensive, be prepared 
for this protocol to take a long time.  The `optimize` commands tries to 
minimize the number of full-length FastRelax simulations it needs to run by 
doing a "preoptimization".  This is just a regular optimization with fewer 
rotamers and (typically) a looser convergence threshold.  The results of this 
simulation, if present, will be used to narrow the range of restraint weights 
sampled during optimization, which will reduce the number of times the 
structure needs to be relaxed.
"""

import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta import *

import numpy as np
import pandas as pd
from math import pi

import sys, os, logging, shutil, stat, json, toml, functools
from pathlib import Path
from appdirs import AppDirs

logger = logging.getLogger('dnp_relax_b')

def main():
    try:
        init_logging()

        import docopt
        args = docopt.docopt(__doc__)

        # Create or load a workspace:
        if args['init']:
            init_workspace(
                    xtal_path=args['<xtal>'],
                    name=args['<workspace>'],
                    overwrite=args['--force'],
                    debug=args['--debug'],
            )
            return

        work = Workspace(args['<workspace>'])

        # Do the simulation/analysis requested by the user:
        if args['optimize']:
            if not args['--skip-preoptimize']:
                run_preoptimize_step(work, overwrite=args['--force'])
            run_optimize_step(work, overwrite=args['--force'])

        if args['relax']:
            init_rosetta(work.config['debug'])
            run_relax_step(work, id=args['<id>'], overwrite=args['--force'])

        if args['finish']:
            pick_best_relaxed_model(work, overwrite=args['--force'])

    except RefusingtoOverwrite as err:
        logger.error(err)
    except KeyboardInterrupt:
        print()

def init_logging():
    # Log everything to stdout.  Individual steps may also setup logging to a 
    # file.
    log = logging.getLogger()
    log.setLevel('DEBUG')

    formatter = logging.Formatter('{asctime}\t{name}\t{levelname}\t{message}', style='{')
    handler = logging.StreamHandler()

    for handler in log.handlers[:]:
        log.removeHandler(handler)

    handler.setFormatter(formatter)
    log.addHandler(handler)

def init_workspace(xtal_path, name='relax_?', overwrite=False, debug=False):
    # Make the root directory:
    xtal_path = Path(xtal_path)
    root = Path(name.replace('?', xtal_path.stem)).resolve()

    if root.exists():
        if overwrite:
            shutil.rmtree(root)
        else:
            raise RefusingtoOverwrite(root)

    root.mkdir()

    # Copy the default configuration file into the workspace:
    with (root/'conf.toml').open('w') as f:
        f.write(f"""\
# See <github.com/toml-lang/toml> for a description of the TOML syntax.
        
# Dramatically reduce the quality of the simulation in order to complete 
# quickly.  In particular, the only one FastRelax iteration is executed and a 
# reduced set of rotamers is used for repacking.
debug = {'true' if debug else 'false'}

# The percent of atoms that should have remain within the average mean squared 
# displacement defined by their B-factors.  This option provides a knob 
# (hopefully one that is meaningful and intuitive) to control how much the 
# model moves during relaxation.  Higher values will produce less movement, as 
# fewer atoms are allowed to move significantly, while lower value will produce 
# more movement.  
preoptimize.percent-within = 90
optimize.percent-within = 90

# How closely the objective (see `percent-within`) needs to be satisfied before 
# the optimization converges.  The defaults specify that 85-95% of the atoms in 
# the model must be consistent with their B-factors.
preoptimize.convergence-tolerance = 10
optimize.convergence-tolerance = 5

# The tightest and loosest restraints that will be considered, respectively.  
# The values refer to the width of a harmonic well that will be applied to each 
# backbone atom, and should be separated by a comma.  You can help the 
# optimization converge faster by specifying a narrower range, but if the 
# actual tightness needed to achieve the given percentage of consistent atoms 
# is outside this range, the simulation will fail and you will need to manually 
# specify a larger range.  Note that `optimize.restraint-range` will be 
# overridden by `optimize.preoptimize-delta` if preoptimization was run.
preoptimize.restraint-range = [0.5, 5.0]
optimize.restraint-range = [0.5, 5.0]

# After the preoptimize step, set the restraint range for the optimize step to 
# the optimal value found during preoptimization plus/minus the specified 
# value.  This overrides `optimize.restraint-range`.
optimize.preoptimize-delta = 0.5
""")

    # Get the paths to everything in the workspace:
    work = Workspace(root)

    # Copy the input structure into the workspace.  I decided to copy rather 
    # than link because copies are more robust and therefore a better record of 
    # what the simulation actually did.  Furthermore, the symlink shouldn't be 
    # necessary to figure out which input was used, because the workspace 
    # should be given a descriptive name.
    shutil.copy(xtal_path, work.unrelaxed_path)

    # Return a recreated the workspace so that the new configuration file will 
    # be properly loaded.
    return Workspace(root)

def run_preoptimize_step(work, overwrite=False):
    step = work.steps['preoptimize']
    step.init(overwrite)

    init_rosetta(test_cycles=True)
    optimize_restraints(
            work.unrelaxed_pose,
            target_percentile=step.config['percent-within'],
            tolerance_percentile=step.config['convergence-tolerance'],
            interval=step.config['restraint-range'],
            extra_rotamers=False,
            status_callback=step.record,
            debug=work.config['debug'],
    )

def run_optimize_step(work, overwrite=False):
    step = work.steps['optimize']
    step.init(overwrite)

    if not work.is_preoptimized:
        interval = step.config['restraint-range']
    else:
        center = work.steps['preoptimize'].optimum['cst_stdev']
        delta = step.config['preoptimize-delta']
        interval = center - delta, center + delta

    init_rosetta(test_cycles=work.config['debug'])
    optimize_restraints(
            work.unrelaxed_pose,
            target_percentile=step.config['percent-within'],
            tolerance_percentile=step.config['convergence-tolerance'],
            interval=interval,
            status_callback=step.record,
            debug=work.config['debug'],
    )

    work.record_optimal_cst_stdev(step.optimum['cst_stdev'])

def run_relax_step(work, id=None, overwrite=False):
    if id is None:
        id = os.environ.get('SLURM_ARRAY_TASK_ID', '0')

    step = work.steps['relax']
    step.init(id, overwrite)

    init_rosetta(test_cycles=work.config['debug'])
    relaxed_pose = core.pose.Pose(work.unrelaxed_pose)
    sfxn = relax_pose(
            relaxed_pose,
            cst_stdev=work.optimal_cst_stdev,
            debug=work.config['debug'],
    )
    relaxed_score = sfxn(relaxed_pose)
    score_diff = relaxed_score - sfxn(work.unrelaxed_pose)

    step.record(
            id,
            relaxed_pose,
            score_reu=relaxed_score,
            score_diff_reu=score_diff,
    )

def pick_best_relaxed_model(work, overwrite=False):
    step = work.steps['relax']

    def score_from_pdb(pdb_path):
        with open(pdb_path) as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith('pose'):
                return float(line.split()[-1])

    scores = {
            Path(x): score_from_pdb(x)
            for x in step.pdb_paths
    }
    if not scores:
        logger.error("no *.pdb files in '{step.root}'")
        return

    models = sorted(scores.keys(), key=lambda x: scores[x])
    best_model = models[0]

    logger.info(f"{len(models)} relaxed models.")
    logger.info(f"Sorted scores:")
    for model in models:
        logger.info(f"{model.name:>10s}: {scores[model]:.3f}")

    logger.info(f"Symlinking 'relaxed.pdb' to '{best_model.relative_to(work.root)}'.")
    work.record_relaxed_pose(best_model, overwrite)


def optimize_restraints(initial_pose, *,
        target_percentile=90, 
        tolerance_percentile=1,
        interval=(0.5, 5.0),
        extra_rotamers=True,
        status_callback=lambda *args: None,
        debug=False,
    ):

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

    def objective(cst_stdev):
        nonlocal iteration
        iteration += 1

        # Relax the pose with the given restraint.
        relaxed_pose = core.pose.Pose(initial_pose)
        sfxn = relax_pose(
                relaxed_pose,
                cst_stdev=cst_stdev,
                extra_rotamers=extra_rotamers,
                debug=debug,
        )

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

        status_callback(
                iteration, 
                cst_stdev,
                relaxed_pose,
                objective=objective,
                score_reu=relaxed_score,
                score_diff_reu=score_diff,
        )

        logger.info(f"Iteration:                        {iteration}")
        logger.info(f"Restraint weight:                 {cst_stdev}")
        logger.info(f"Mean squared distance (actual):   {df.dist2_actual.mean()}")
        logger.info(f"Mean squared distance (B-factor): {df.dist2_b_factor.mean()}")
        logger.info(f"Objective:                        {objective}")
        logger.info(f"Score (REU):                      {relaxed_score}")
        logger.info(f"Score improvement (REU):          {score_diff}")

        return objective

    logger.info(f"Target percent within:            {target_percentile}%")
    logger.info(f"Convergence tolerance:            {tolerance_percentile}%")
    logger.info(f"Min restraint weight:             {interval[0]}")
    logger.info(f"Max restraint weight:             {interval[1]}")

    x, results = brentq(
            objective,
            interval[0], interval[1],
            xtol=tolerance_percentile / 100,
            full_output=True,
    )

    logger.info(f"Converged:                        {results.converged}")
    logger.info(f"Best restraint weight:            {results.root}")
    logger.info(f"Iterations:                       {results.iterations}")
    logger.info(f"Function calls:                   {results.function_calls}")

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


def init_rosetta(test_cycles=False, constant_seed=False):
    flags = [
            '-relax:dna_move on',
            '-relax:coord_cst_stdev 0.5',  # default: 0.5
            '-relax:constrain_relax_to_start_coords on',
            '-relax:ramp_constraints off',
            '-run:test_cycles', 'on' if test_cycles else 'off',
            '-dna:specificity:exclude_dna_dna off',
            '-out:levels core.pack.rotamer_set.RotamerSet_.extra_rotamers:500',
    ]

    if constant_seed:
        flags += [
            '-constant_seed',
        ]

    pyrosetta.init(' '.join(flags), set_logging_handler='logging')

def relax_pose(pose, cst_stdev=1.0, extra_rotamers=True, debug=False):
    # The tightness of the coordinate restraints can only be set via the 
    # command-line.  Welcome to Rosetta...
    basic.options.set_real_option('relax:coord_cst_stdev', cst_stdev)

    sfxn = load_score_function()
    tf = load_task_factory(extra_rotamers, debug)
    mm = load_move_map()
    relax = load_fast_relax(sfxn, tf, mm)

    relax.apply(pose)

    return sfxn

def load_score_function():
    return core.scoring.ScoreFunctionFactory.create_score_function('ref2015')

def load_task_factory(extra_rotamers=True, debug=True):
    tf = core.pack.task.TaskFactory()

    # FastRelax automatically adds the `IncludeCurrent` task operation, so 
    # adding it here is redundant, but hopefully helpful in terms of clarity.
    tf.push_back(core.pack.task.operation.IncludeCurrent())
    tf.push_back(core.pack.task.operation.RestrictToRepacking())

    if extra_rotamers and not debug:
        load_extra_rotamers(tf)

    return tf

def load_extra_rotamers(tf):
    from pyrosetta.rosetta.core.pack.task.operation import (
            OperateOnResidueSubset,
    )
    from pyrosetta.rosetta.core.select.residue_selector import (
            ResiduePropertySelector,
            AndResidueSelector as And,
            NotResidueSelector as Not,
    )

    protein_sele = ResiduePropertySelector(core.chemical.PROTEIN)
    dna_sele = ResiduePropertySelector(core.chemical.DNA)
    interface_sele = load_interface_sele(protein_sele, dna_sele)
    noninterface_sele = Not(interface_sele)

    task_op_args = [
            (protein_sele, interface_sele,    load_protein_interface_rlt),
            (protein_sele, noninterface_sele, load_protein_noninterface_rlt),
            (dna_sele,     interface_sele,    load_dna_interface_rlt),
            (dna_sele,     noninterface_sele, load_dna_noninterface_rlt),
    ]

    for chain_sele, interface_sele, rlt in task_op_args:
        sele = And(chain_sele, interface_sele)
        tf.push_back(OperateOnResidueSubset(rlt(), sele))

def load_interface_sele(protein_sele, dna_sele):
    return core.select.residue_selector.\
            InterGroupInterfaceByVectorSelector(protein_sele, dna_sele)

def load_protein_interface_rlt():
    ex = core.pack.task.operation.ExtraRotamersGenericRLT()

    # Use extra rotamers for all χ angles near the interface: Arg and Lys are 
    # important for DNA interfaces, so we want to sample those sidechains 
    # finely.
    ex.ex1(True)
    ex.ex2(True)
    ex.ex3(True)
    ex.ex4(True)
    ex.extrachi_cutoff(1)

    return ex

def load_protein_noninterface_rlt():
    ex = core.pack.task.operation.ExtraRotamersGenericRLT()
    ex.ex1(True)
    ex.ex2(True)
    ex.extrachi_cutoff(1)
    return ex

def load_dna_interface_rlt():
    ex = core.pack.task.operation.ExtraRotamersGenericRLT()
    ex.exdna_sample_level(core.pack.task.NO_EXTRA_CHI_SAMPLES) # 1 rot/nt
    return ex

def load_dna_noninterface_rlt():
    return load_dna_interface_rlt()
    
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

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.config_path = self.root / 'conf.toml'
        self.relaxed_path = self.root / 'relaxed.pdb'
        self.unrelaxed_path = self.root / 'unrelaxed.pdb'
        self.cst_stdev_path = self.root / 'cst_stdev'

        # Create representations for the various steps of the protocol.  This 
        # has to be the last thing done, because these step will set themselves 
        # up assuming that the workspace is ready to go.
        self.steps = {
                'preoptimize': OptimizeStep(self, 'preoptimize'),
                'optimize': OptimizeStep(self, 'optimize'),
                'relax': RelaxStep(self),
        }

    @property
    @functools.lru_cache(None)
    def config(self):
        return toml.load(self.config_path)

    @property
    def unrelaxed_pose(self):
        return pose_from_pdb(str(self.unrelaxed_path))

    def record_relaxed_pose(self, path, overwrite=False):
        if self.relaxed_path.exists():
            if overwrite:
                self.relaxed_path.unlink()
            else:
                raise RefusingtoOverwrite(self.relaxed_path)

        self.relaxed_path.symlink_to(
                Path(path).relative_to(self.root))

    def record_optimal_cst_stdev(self, cst_stdev):
        with self.cst_stdev_path.open('w') as f:
            f.write(f'{cst_stdev}\n')

    @property
    def optimal_cst_stdev(self):
        if not self.cst_stdev_path.exists():
            raise NotOptimized(self.cst_stdev_path)

        with self.cst_stdev_path.open() as f:
            return float(f.read())

    @property
    def is_preoptimized(self):
        return self.steps['preoptimize'].trajectory_path.exists()


class Step:

    def __init__(self, work, name):
        self.work = work
        self.name = name
        self.root = self.work.root / self.name

    @property
    @functools.lru_cache(None)
    def config(self):
        return self.work.config.get(self.name, {})


class OptimizeStep(Step):

    def __init__(self, work, name):
        super().__init__(work, name)

        # Track each step in the optimization process.
        self.trajectory_path = self.root / 'trajectory.json'

        if not self.trajectory_path.exists():
            self.trajectory = {}
        else:
            with self.trajectory_path.open() as f:
                self.trajectory = json.load(f)

    def init(self, overwrite=False):
        # Make a directory for this step.
        if self.root.exists():
            if overwrite:
                shutil.rmtree(self.root)
            else:
                raise RefusingtoOverwrite(self.root)

        self.root.mkdir()

        # Log results from this simulation to a file.
        log = logging.getLogger()
        handler = logging.FileHandler(self.root / 'log')
        log.addHandler(handler)

    def record(self, iteration, x, pose, **kwargs):
        # Write the PDB.
        path = self.root / f'iteration_{iteration}.pdb'
        pose.dump_pdb(str(path))

        # Add to the JSON log.
        if iteration in self.trajectory:
            logger.error(f"overwriting trajectory information for iteration {iteration}: {trajectory[iteration]}")

        self.trajectory[iteration] = {
                'path': str(path),
                'cst_stdev': x,
                **kwargs,
        }

        with self.trajectory_path.open('w') as f:
            json.dump(self.trajectory, f)

    @property
    def optimum(self):
        frames = sorted(self.trajectory.values(), key=lambda x: abs(x['objective']))
        return frames[0]

class RelaxStep(Step):

    def __init__(self, work):
        super().__init__(work, 'relax')

        self.scores_path = self.root / 'scores.json'

        if not self.scores_path.exists():
            self.scores = {}
        else:
            with self.scores_path.open() as f:
                self.scores = json.load(f)

    def init(self, id, overwrite=False):
        # Make a directory for this step.
        if self.root.exists() and overwrite:
            shutil.rmtree(self.root)

        self.root.mkdir(exist_ok=True)

        # Log results from this simulation to a file.
        log = logging.getLogger()
        formatter = logging.Formatter('{asctime}\t{name}\t{levelname}\t{message}', style='{')
        handler = logging.FileHandler(self.root / f'{id}.log')

        handler.setFormatter(formatter)
        log.addHandler(handler)

    def record(self, id, pose, **kwargs):
        id = str(id)
        pdb_path = self.root / f'{id}.pdb'
        json_path = self.root / f'{id}.json'

        # Write the PDB.
        if pdb_path.exists():
            logger.warning(f"overwriting '{pdb_path}'")
        logger.info(f"writing relaxed pose to '{pdb_path}'")

        pose.dump_pdb(str(pdb_path))

        # Write any keyword arguments to a JSON file.
        if json_path.exists():
            logger.warning(f"overwriting '{json_path}'")
        logger.info(f"writing extra information to '{json_path}'")

        with json_path.open('w') as f:
            x = {'path': pdb_path, 'id': id, **kwargs}
            json.dump(x, f)

    @property
    def pdb_paths(self):
        yield from self.root.glob('*.pdb')



class RefusingtoOverwrite(Exception):

    def __init__(self, path):
        super().__init__(f"'{path}' already exists.  Use `-f` to overwrite it.")

class NotOptimized(Exception):

    def __init__(self, path):
        super().__init__(f"no optimization results found: '{path}'")
