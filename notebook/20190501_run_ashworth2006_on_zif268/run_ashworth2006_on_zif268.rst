*****************************
Run [Ashworth2006]_ on Zif268
*****************************

My goal is to run as many of the existing protein/DNA interface design 
algorithms on as many different scaffolds as possible, but I want to start with 
the simplest algorithm and the most canonical scaffold.  [Ashworth2006]_ is 
basically just packing and minimizing, and Zif268 is the scaffold I'm using to 
validate my high-throughput screen.

Note that [Thyme2014]_ describes how to run this protocol using RosettaScripts, 
but I'll have to adapt those instructions because I definitely want to use 
PyRosetta.  I like the idea of treating Rosetta like a true library, and 
potentially being able to build nice frontends for particular design projects 
using python.

Comments
========

``INCLUDE_DNA_DNA``
-------------------
From past experience, I know that the ``INCLUDE_DNA_DNA`` flag is required when 
allowing DNA to pack.  However, I had to read the code to figure out exactly 
what it does:

- The ``INCLUDE_DNA_DNA`` flag sets the ``exclude_DNA_DNA`` energy method 
  option, so most of the references in the code base are to 
  ``exclude_DNA_DNA``.  
  
- Directly setting the ``exclude_DNA_DNA`` energy method option in PyRosetta 
  does not seem to have any effect on scoring, regardless of whether I recreate 
  the score function or clear the energies cache.  However, both putting the 
  ``INCLUDE_DNA_DNA`` flag in the weights file and adding 
  ``-dna:specificity:exclude_dna_dna off`` to the command line do affect 
  scoring.  I'm not sure why, but it seems like the initial setting is what 
  matters.

  I'd like to not burden the user with this detail, so adding the command-line 
  argument is probably the way to go.
  
- The score terms that seems to be affected by ``exclude_DNA_DNA`` (based on 
  ack-ing ``src/core/scoring`` for energy methods and looking at the terms 
  declared by the corresponding creators) are:

   - ``fa_atr``
   - ``fa_rep``
   - ``fa_sol``

   - ``fa_elec``
	- ``fa_elec_bb_bb``
	- ``fa_elec_bb_sc``
	- ``fa_elec_sc_sc``
	- ``fa_intra_elec``
   - ``fa_grpelec``

   - ``hbond_lr_bb``
	- ``hbond_sr_bb``
	- ``hbond_bb_sc``
	- ``hbond_sr_bb_sc``
	- ``hbond_lr_bb_sc``
	- ``hbond_sc``
   - ``hbond_wat``
   - ``wat_entropy``
   - ``hbond_intra``
	- ``hbond``

   - ``facts_elec``
   - ``facts_solv``
   - ``facts_sasa``
   - ``gb_elec``
   - ``fa_vdw_tinker``
   - ``multipole_elec``
   - ``fa_sasa``
   - ``buried_unsat_penalty``

- In most cases, ``exclude_DNA_DNA`` seems to just turn off terms for DNA/DNA 
  interactions, if it's not disabled.

- When scoring 1BNA without ``INCLUDE_DNA_DNA``, only the ``lk_ball_wtd`` term 
  has a non-zero value.

- There is also an ``INCLUDE_HB_DNA_DNA`` option, which corresponds to 
  ``hbond_options()->exclude_DNA_DNA()``.  Presumably this causes the H-bonding 
  term to be calculated between DNA nucleotides.  This option is automatically 
  set by ``INCLUDE_DNA_DNA``, so it should only be necessary to set if you want 
  to score intra-DNA H-bonds, but nothing else.  
  
  Actually, attempting to use either one of the ``DNA_DNA`` flags without the 
  other triggers an assertion error.  So these are clearly meant to be the same 
  thing, and just are not very DRY for whatever reason.

Packing
=======
- In order to repack DNA, the ``exdna_sample_level`` option needs to be set.  
  This option controls how many rotamers are generated per base.  It's a bit 
  confusing though, because although the option uses the "extra chi" enum, it's 
  behavior has nothing to do with standard deviations or half-steps:

   ===============================  =====  ===============
   Enum                             Value      Rotamers/nt
   ===============================  =====  ===============
   ``NO_EXTRA_CHI_SAMPLES``             0                1
   ``EX_ONE_STDDEV``                    1               10
   ``EX_ONE_HALF_STEP_STDDEV``          2               25
   ``EX_TWO_FULL_STEP_STDDEVS``         3               50
   ``EX_TWO_HALF_STEP_STDDEVS``         4              100
   ``EX_FOUR_HALF_STEP_STDDEVS``        5              250
   ``EX_THREE_THIRD_STEP_STDDEVS``      6              500
   ``EX_SIX_QUARTER_STEP_STDDEVS``      7  ~20 (see below)
   ===============================  =====  ===============

- ``EX_SIX_QUARTER_STEP_STDDEVS`` is a special case:  It causes rotamers to be 
  loaded from a file in the database: ``rotamer/dna/VQ-DNA-64.rotlib``.  Note 
  that there are other DNA rotamer libraries in that directory, but as far as I 
  can tell only the "64" one can be used (without altering the code).  

- The rest of the "ex" settings, e.g. ``-ex1 -ex2``, don't affect DNA.

- At least for 1BNA, repacking alone is not enough to get reasonable scores.  
  Even when packing with the maximal number of rotamers, the final score is ~81 
  REU worse than the starting structure, mostly due to repulsive interactions 
  in the base-pair H-bonding interface.

Minimization
============
- Minimization seems to work just fine with DNA.

- The caveat is to be careful about minimizing after repacking.  Repacking can 
  introduce pretty significant clashes, more than it seems to with proteins, in 
  which case subsequent minimization can explode the structure.  Fixed backbone 
  minimization after packing helps, but does not solve the problem.  Coordinate 
  restraints or a soft repulsive term would surely help as well.

Relax
=====
I'm trying to find a good way to relax DNA structures.

- The PDB_redo version of 1BNA scores significantly better than the version in 
  the PDB, despite moving only very slightly.  So clearly, it should be 
  possible to get score improvements with minimal change in structure.

- I can't get the same behavior out of minimization.  If I don't let the 
  backbone minimize, the score can hardly improve at all.  If I let the 
  backbone minimize, it moves way too much lot.  Cartesian minimization keeps 
  things more in place, but it's slow and it still moves things enough to 
  really change the structure.

- ``-relax:constrain_relax_to_start_coords on`` seems to work, while 
  ``relax.constrain_relax_to_start_coords()`` seems to have no effect.  This 
  code is so fucking bad...

   - Ok, in python, you need to explicitly set ``relax.constrain_coords()`` in 
     addition to the setting to constrain to either starting of native coords.  
     So fucking bad...

- The ``relax.ramp_down_constraints()`` method doesn't really seem to do 
  anything.  However, both the ``-relax:ramp_constraints off`` flag and fast 
  relax protocols that don't ramp constraints (the third number on the 
  ``ramp_repack_min`` lines) do prevent constraint ramping.  I looked at the 
  code a bit, but couldn't see why this would be the case.  For now, just use 
  one of the methods that works.
  
  Maybe.  If I change the fast relax protocol to have a constraint weight of 1 
  the whole time, the constraints are respected.  If I use the standard 
  protocol with ``relax.ramp_down_constraints(False)``, the backbone moves a 
  lot and I ends up with an extremely high score, due to the coordinate term.

- Doesn't seem necessary to restrain sidechain atoms.  I think the base pairs
  really just don't have much opportunity to move, if the backbone is 
  restrained.

- Relaxing with 50 rotamers/nt didn't perform any better than relaxing with 25 
  rotamers/nt

- The IncludeCurrent task operation doesn't seem to do anything.  In 
  particular, it doesn't seem to affect the number of rotamers considered at 
  each packing step.  Since I found that including the current rotamer was 
  pretty important with just packing, I think relax is just automatically 
  including the current rotamers.  It might be worth trying to read the code to 
  confirm this.

  - Yes, confirmed.  See ``FastRelax::apply()``.  The ``IncludeCurrent`` task 
    operation is unconditionally added.

  - So looking at the scores below, it seems that repacking (at least for this 
    structure) is not really doing anything: the crystal structure rotamer is 
    the best.  I'm sure this is because all the other rotamers clash horribly.

- The looser the coordinate constraints are: the more the backbone moves, and 
  the better the score gets.  I'm not really interested in sampling DNA 
  backbone flexibility during relaxation, though.  So I think what I need to 
  find is how tight to make the restraints such that I get backbone movement 
  that's on par with the average crystal structure error.

  Right now I'm doing this by comparing relaxed structures to 1BNA and 1BNA 
  redo.  Specifically I look at the backbone sugar, and judge if the relaxed 
  atoms are in line with those from the two crystal structures.  This could be 
  more rigorous: I could scan different constraint sigmas to find the one that 
  gives me the same backbone RMSD as the average coordinate error (I think 
  ~0.3Å).

  One question is whether it would be better to use torsion or cartesian 
  minimization.  Cartesian minimization is better at not blowing up the 
  structure in general, as so give smaller backbone perturbations when given 
  more freedom to move.  This does not appear to be the case, however.  At 
  σ=5.0, both torsion and cartesian minimization give similar amounts of 
  movement, and produce structures that are much more similar to each other 
  than to the crystal structures.  

  It seems that the strength of the coordinate restraints is really the most 
  important parameter. 

=========  ======  =======  =====  ===========  =========  ====================
Structure  Relax?  Rot/nt?  Cart?  Restraints?      Score  Notes
=========  ======  =======  =====  ===========  =========  ====================
1BNA                                            -12.94055
1BNA         X                                  -32.16100
1BNA redo                                       -27.46022
1BNA redo    X                                  -41.68700
---------  ------  -------  -----  -----------  ---------  --------------------
1BNA         X           1                      -32.16100
1BNA         X          10                      -32.16100
1BNA         X          25                      -32.16100
1BNA         X          50                      -31.88394 
---------  ------  -------  -----  -----------  ---------  --------------------
1BNA         X           1    X                 -31.72935
1BNA         X          10    X                 -31.72965
1BNA         X          25    X                 -33.05246
1BNA         X          50    X                 -31.72839
---------  ------  -------  -----  -----------  ---------  --------------------
1BNA         X           1                 0.5  -32.16100
1BNA         X           1    X            0.5  -38.60980
1BNA         X           1                 1.0  -46.05915
1BNA         X           1    X            1.0  -52.62182  Not bad
1BNA         X           1                 5.0  -68.79233  Backbone clearly translated.
1BNA         X           1    X            5.0  -84.23563  Backbone clearly translated.
=========  ======  =======  =====  ===========  =========  ====================


References
==========
.. [Ashworth2006] :doi:`10.1038/nature04818`
.. [Thyme2014] :doi:`10.1007/978-1-62703-968-0_17`

