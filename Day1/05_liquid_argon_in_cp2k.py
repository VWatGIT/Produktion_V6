#!/usr/bin/python

""".. _liquid_argon_in_cp2k:

Exercise 1.4: Liquid argon simulation using ASE and CP2K
=========================================================

In this tutorial, we learn how to simulate liquid argon at 85 K using
Born-Oppenheimer ab-initio molecular dynamics. We use ASE (Atomic Simulation
Environment) to set up the system and run the MD, while CP2K computes
the DFT forces at each timestep.

First, we import our Python modules.
"""

import ase.io
from ase.calculators.cp2k import CP2K
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase import units

# %%
# Loading the initial configuration
# ---------------------------------
#
# We load the starting configuration of 108 argon atoms from an XYZ file,
# set the cubic simulation box, and apply periodic boundary conditions.

atoms = ase.io.read("system.xyz")
atoms.set_cell([17.0742, 17.0742, 17.0742])
atoms.pbc = True

# %%
# Setting up the DFT force calculation
# ------------------------------------
#
# Next, we configure CP2K as our DFT calculator. We only need to specify
# the DFT settings here -- ASE will handle the molecular dynamics.
#
# We pass additional CP2K settings via the ``inp`` parameter. These control
# the SCF convergence and orbital transformation settings.

inp = """
&FORCE_EVAL
    &DFT
        &SCF
            &OT
                MINIMIZER DIIS
                PRECONDITIONER FULL_ALL
            &END OT
        &END SCF
    &END DFT
&END FORCE_EVAL
"""

# %%
# Now we create the CP2K calculator and attach it to our atoms.
# We specify the basis set, pseudopotential, and exchange-correlation functional.
#
# .. admonition:: About basis sets
#    :class: info
#
#    A basis set is a set of functions used to represent the electronic
#    wave function.
#
# .. admonition:: About pseudopotentials
#    :class: info
#
#    The pseudopotential approximates the nucleus-valence electron and
#    core electron-valence electron interaction.
#
# .. admonition:: About exchange correlation functionals
#    :class: info
#
#    The exchange correlation functional approximates the electronic
#    exchange and correlation energy from the electron density. Here we
#    use the Perdew-Burke-Ernzerhof functional (PBE).
#
# Make sure you have the files ``BASIS_MOLOPT`` and ``GTH_POTENTIALS``
# in your working directory (or accessible to CP2K).

atoms.calc = CP2K(
    inp=inp,
    basis_set="DZVP-MOLOPT-SR-GTH",
    pseudo_potential="GTH-PBE-q8",
    potential_file="GTH_POTENTIALS",
    xc="PBE",
    cutoff=400 * units.Rydberg,
    command="cp2k.psmp -s",
)

# %%
# Initializing velocities
# -----------------------
#
# Before starting the MD, we need to assign initial velocities to all atoms
# according to a Maxwell-Boltzmann distribution at the target temperature.

MaxwellBoltzmannDistribution(atoms, temperature_K=85)

# %%
# Setting up the molecular dynamics
# ---------------------------------
#
# We use Langevin dynamics, which couples the system to a heat bath at
# constant temperature (NVT ensemble). The ``friction`` parameter controls
# how strongly the thermostat acts on the atoms.
#
# .. admonition:: About the thermostat
#    :class: info
#
#    A thermostat modifies particle velocities during the simulation to
#    keep the temperature constant. Langevin dynamics adds a friction term
#    and random forces to the equations of motion, effectively coupling
#    the system to a heat bath.

dyn = Langevin(atoms, timestep=10 * units.fs, temperature_K=85, friction=0.01 / units.fs)

# %%
# Writing output
# --------------
#
# We attach a trajectory writer that saves the atomic positions, velocities,
# and forces at every timestep. This is the ASE equivalent of the CP2K
# PRINT section for trajectories, velocities, and forces.

traj = Trajectory("Argon_Simulation.traj", "w", atoms)
dyn.attach(traj.write, interval=1)

# %%
# Running the simulation
# ----------------------
#
# Now we run 5000 MD steps. At each step, ASE calls CP2K to compute
# the DFT forces, then advances the atomic positions using the
# Langevin integrator.
#
# .. note::
#
#    This simulation generates the trajectory file ``Argon_Simulation.traj``
#    that will be analysed at the start of Day 2 (:ref:`analyse_cp2k_trajectory`).

dyn.run(5000)

# %%
# Estimating the runtime
# ----------------------
#
# Each MD step requires a full DFT self-consistent field (SCF) cycle for
# 108 atoms — this is orders of magnitude more expensive than a classical
# force-field evaluation.
#
# .. warning::
#
#    Running 5000 DFT steps will take **several hours** on a typical workstation.
#    Start the simulation now and let it run in the background.
#
# CP2K writes its full output to a file called ``cp2k.out`` in your working
# directory.
#
# .. admonition:: Task 1
#
#    After the first few steps have completed, open ``cp2k.out`` and find the
#    wall-time reported at the end of a single MD step.  Look for a line like::
#
#    Multiply this wall time per step by 5000 to estimate the total runtime.
