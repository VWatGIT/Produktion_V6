#!/usr/bin/python

""".. _analyse_cp2k_trajectory:

Exercise 2.0: Analysing the liquid argon CP2K trajectory
=========================================================

Welcome to Day 2.  Before training a machine learning potential we first
inspect the DFT trajectory generated in Exercise 1.4
(:ref:`liquid_argon_in_cp2k`).  Checking equilibration, temperature
stability, and density gives us confidence that the data is physically
reasonable and ready for training.

First, we import our Python modules.
"""

import matplotlib.pyplot as plt
import ase.io

# %%
# Loading the trajectory
# ----------------------
#
# We read back the ASE trajectory written by the CP2K MD simulation.
# Each frame already carries the DFT potential energy and forces
# computed by CP2K.

trajectory = ase.io.read("Argon_Simulation.traj", index=":")
print(f"Loaded {len(trajectory)} frames")

# %%
# Potential energy over time
# --------------------------
#
# .. admonition:: Task 1
#
#    Plot the **potential energy** as a function of MD step.
#    Use ``frame.get_potential_energy()`` to extract the energy from each
#    frame.
#
#    Has the system equilibrated?  If not, how many frames should be
#    discarded before training?

pass  # TODO: extract and plot the potential energy for each frame

# %%
# Temperature over time
# ---------------------
#
# .. admonition:: Task 2
#
#    Plot the **instantaneous temperature** as a function of MD step.
#    Use ``frame.get_temperature()`` (which computes temperature from
#    the kinetic energy of the frame).
#
#    Does the Langevin thermostat keep the temperature around 85 K?

pass  # TODO: extract and plot the instantaneous temperature for each frame
