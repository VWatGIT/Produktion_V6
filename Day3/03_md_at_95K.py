#!/usr/bin/python

""".. _md_at_95K:

Exercise 3.3: MD of argon at 95 K — testing the limits
========================================================

In this tutorial we test how well our machine learning potential
generalises beyond its training conditions. The Apax model was trained
on DFT data at 85 K. We now run MD at 95 K and compare the results.

This exercise illustrates an important limitation of ML potentials:
they are most reliable within (or close to) the thermodynamic conditions
represented in the training data.

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
# We load the same starting configuration and attach the Apax calculator
# trained at 85 K.

PROJECT_PATH = Path("../../../solutions")

atoms = ase.io.read(PROJECT_PATH / "dft/liquid_argon_85K/system.xyz")
atoms.set_cell([17.0742, 17.0742, 17.0742])
atoms.pbc = True

model_dir = PROJECT_PATH / "apax" / "Argon_Simulation"
calc = ASECalculator(model_dir=model_dir)
atoms.calc = calc

# %%
# Running MD at 95 K with the 85 K model
# ----------------------------------------
#
# We set the target temperature to 95 K instead of 85 K.
#
# .. admonition:: Task 1
#
#    Run an MD simulation at 95 K using the model trained at 85 K.
#    Use the same simulation parameters as in Exercise 3.1 but change
#    the target temperature to 95 K. Save the trajectory to
#    ``Argon_ML_95K_model85K.traj``.

# TODO: set up and run MD at 95 K, save trajectory

# %%
# Analysing the results
# ----------------------
#
# .. admonition:: Task 2
#
#    Read the 95 K trajectory and plot the potential energy and
#    temperature over time. Compare them to your 85 K results from
#    Exercise 3.1. Do you notice any differences or instabilities?

# TODO: read trajectory and plot energy and temperature

# %%
# Comparing RDFs
# ---------------
#
# .. admonition:: Task 3
#
#    Compute the radial distribution function for the 95 K trajectory
#    and overlay it with the RDF from the 85 K simulation. How do the
#    peak positions and heights compare?
