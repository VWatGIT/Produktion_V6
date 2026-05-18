#!/usr/bin/python

""".. _preprocess_dft:

Exercise 2.1: Preprocessing DFT Data
======================================

In this tutorial we read the DFT trajectory from the previous day and
split it into training, validation, and test sets for Apax.

The trajectory file ``Argon_Simulation.traj`` was produced by the Day 1
liquid argon exercise (:ref:`liquid_argon_in_cp2k`). It is an ASE
trajectory that already contains positions, energies, and forces in
``eV`` and ``eV/Å``, so no unit conversion is needed.

We start with importing our Python modules.
"""

from pathlib import Path
import numpy as np
from ase.io import read, write

# %%
# Then we set our project path. Replace this with your own project path.

PROJECT_PATH = Path("../../../solutions")

# %%
# Reading the ASE trajectory
# ---------------------------
#
# We read the trajectory file that was written by ASE during the Day 1
# CP2K MD simulation. Each frame already carries the DFT energy and
# forces computed by CP2K.

trajectory = read(
    PROJECT_PATH / "dft/liquid_argon_85K/Argon_Simulation.traj",
    index=":",
)

print(f"Loaded {len(trajectory)} frames")
print(f"First frame: energy = {trajectory[0].get_potential_energy():.6f} eV")
print(f"             forces shape = {trajectory[0].get_forces().shape}")

# %%
# Converting to extxyz
# ---------------------
#
# The ASE ``.traj`` format stores results attached to a calculator. To
# write standalone files that Apax can read, we copy the energy and
# forces into the ``info`` and ``arrays`` dictionaries of each frame and
# write them out in the ``extxyz`` format.


def traj_to_extxyz(frames):
    """Copy calculator results into info/arrays for extxyz export."""
    out = []
    for frame in frames:
        atoms = frame.copy()
        atoms.info["energy"] = frame.get_potential_energy()
        atoms.arrays["forces"] = np.array(frame.get_forces())
        out.append(atoms)
    return out


trajectory_extxyz = traj_to_extxyz(trajectory)

# %%
# Splitting the trajectory
# -------------------------
#
# Before training a machine learning potential, we divide our trajectory
# into three non-overlapping subsets. We skip the first 500 frames to
# discard the initial equilibration phase, then take 1000 frames each for
# **training**, **validation**, and **testing**.
#
# .. admonition:: Task 1
#
#    Split the trajectory into training, validation, and test sets.
#    Skip the first 500 frames (equilibration), then take 1000 frames
#    for each set. Write them to ``train.xyz``, ``validate.xyz``, and
#    ``test.xyz`` in extxyz format.

out_dir = PROJECT_PATH / "apax"
out_dir.mkdir(exist_ok=True)
# TODO: split the trajectory and write train.xyz, validate.xyz, test.xyz

# %%
# Sanity check
# ------------
#
# Reload the saved files and verify they contain the expected data.

check = read(out_dir / "train.xyz", index=0, format="extxyz")
print(f"has energy: {'energy' in check.info}")
print(f"has forces: {'forces' in check.arrays}")
print(f"energy = {check.info.get('energy', 'MISSING')}")
