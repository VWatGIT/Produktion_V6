#!/usr/bin/python

""".. _md_with_apax:

Exercise 3.1: Molecular dynamics with Apax
============================================

In this tutorial we use the Apax model trained on Day 2 to run a
molecular dynamics simulation of liquid argon at 85 K. The setup is
almost identical to the DFT MD from Day 1 — we simply swap the CP2K
calculator for the Apax ASE calculator.

First, we import our Python modules.
"""

import matplotlib.pyplot as plt
import ase.io
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase import units
from apax.md import ASECalculator
from pathlib import Path

# %%
# Setting up the system
# ----------------------
#
# We load the same starting configuration of 108 argon atoms that we
# used on Day 1.

PROJECT_PATH = Path("../../../solutions")

atoms = ase.io.read(PROJECT_PATH / "dft/liquid_argon_85K/system.xyz")
atoms.set_cell([17.0742, 17.0742, 17.0742])
atoms.pbc = True

# %%
# Loading the ML potential
# -------------------------
#
# Instead of CP2K we now use the Apax model trained on Day 2 as the
# force calculator. The interface is the same ASE calculator API, so
# the rest of the simulation code does not change.

model_dir = PROJECT_PATH / "apax" / "Argon_Simulation"
calc = ASECalculator(model_dir=model_dir)
atoms.calc = calc

# %%
# Let us verify the calculator works by computing a single energy.

E = atoms.get_potential_energy()
print(f"Initial potential energy: {E:.4f} eV")

# %%
# Initializing velocities
# -------------------------
#
# We assign initial velocities according to a Maxwell-Boltzmann
# distribution at 85 K.

MaxwellBoltzmannDistribution(atoms, temperature_K=85)

# %%
# Setting up molecular dynamics
# ------------------------------
#
# We use Langevin dynamics (NVT ensemble) just as on Day 1.
#
# .. admonition:: Task 1
#
#    Set up a Langevin integrator with a timestep of 10 fs, a target
#    temperature of 85 K, and a friction coefficient of
#    0.01 fs\ :sup:`-1`. Attach a trajectory writer and run 5000 steps.
#    Compare the wall-clock time to the DFT MD from Day 1.

# TODO: set up Langevin dynamics, attach trajectory writer, and run

# %%
# Analysing the trajectory
# -------------------------
#
# After the simulation completes, we read back the trajectory.
#
# .. admonition:: Task 2
#
#    Read the trajectory file ``Argon_ML_85K.traj`` and plot the
#    **potential energy** as a function of MD step. Compare the result
#    to the DFT MD from Day 1 — is the energy in the same range?

# TODO: read the trajectory and plot the potential energy vs MD step

# %%
#
# .. admonition:: Task 3
#
#    Plot the **instantaneous temperature** as a function of MD step.
#    Does the thermostat keep the temperature around 85 K?

# TODO: plot the instantaneous temperature vs MD step

# %%
#
# .. admonition:: Task 4
#
#    Compute the **density** of the system for each frame and plot it
#    over time. Compare to the experimental liquid argon density at
#    85 K (~1400 kg/m\ :sup:`3`).

# TODO: compute and plot the density over time
