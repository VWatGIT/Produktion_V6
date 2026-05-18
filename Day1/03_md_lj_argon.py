#!/usr/bin/python

""".. _md_lj_argon:

Exercise 1.2: Molecular dynamics of liquid argon with ASE
==========================================================

In this tutorial we run a molecular dynamics (MD) simulation of liquid argon
at 85 K using the Lennard-Jones potential we fitted in the previous tutorial.
We use ASE (Atomic Simulation Environment) to set up the system, compute
forces, and integrate the equations of motion.

First, we import our Python modules.
"""

import numpy as np
import matplotlib.pyplot as plt

import ase.io
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase import units

# %%
# Loading the initial configuration
# ----------------------------------
#
# We load the starting configuration of 108 argon atoms from an XYZ file.
# This is the same configuration we will use in the CP2K tutorial
# later today (:ref:`2_body_potential`).
# We set the cubic simulation box and apply periodic boundary conditions.

atoms = ase.io.read("system.xyz")
atoms.set_cell([17.0742, 17.0742, 17.0742])
atoms.pbc = True

print(f"Number of atoms: {len(atoms)}")
print(f"Cell: {atoms.cell.lengths()}")

# %%
# The equations of motion
# ------------------------
#
# In molecular dynamics we solve Newton's second law for each atom *i*:
#
# .. math::
#
#    m_i \, \mathbf{a}_i = \mathbf{F}_i
#    = -\nabla_i \sum_{j > i} V(r_{ij})
#
# The force on each atom is the negative gradient of the total potential
# energy, which is a sum over all pair interactions.
#
# The positions and velocities are propagated in time using the
# **Velocity Verlet** algorithm:
#
# .. math::
#
#    \mathbf{v}(t + \tfrac{\Delta t}{2}) &= \mathbf{v}(t)
#    + \tfrac{\Delta t}{2}\,\frac{\mathbf{F}(t)}{m} \\
#    \mathbf{r}(t + \Delta t) &= \mathbf{r}(t)
#    + \Delta t\;\mathbf{v}(t + \tfrac{\Delta t}{2}) \\
#    \mathbf{v}(t + \Delta t) &= \mathbf{v}(t + \tfrac{\Delta t}{2})
#    + \tfrac{\Delta t}{2}\,\frac{\mathbf{F}(t + \Delta t)}{m}
#
# This algorithm is time-reversible, symplectic, and conserves energy
# well over long simulations.
#
# Setting up the Lennard-Jones calculator
# ----------------------------------------
#
# We use ASE's built-in ``LennardJones`` calculator.  The parameters
# are those we fitted to experimental data in the previous tutorial.
# ASE expects *epsilon* in eV and *sigma* in Angstrom.
#
# .. seealso::
#
#    - `ASE calculators overview <https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html>`_
#    - `ASE LennardJones calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/lj.html>`_
#    - `Attaching a calculator to an Atoms object <https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.calc>`_
#
# .. admonition:: Task 1
#
#    Set up the LJ calculator for argon with the following parameters:
#
#    - :math:`\varepsilon/k_B = 120\,\text{K}` (convert to eV using
#      :math:`k_B = 8.617333 \times 10^{-5}\,\text{eV/K}`)
#    - :math:`\sigma = 3.4\,\text{\AA}`
#    - cutoff at :math:`2.5\,\sigma`
#
#    Assign the calculator to the ``atoms`` object using ``atoms.calc = ...``.
#
#    .. hint::
#
#       The ``LennardJones`` constructor accepts keyword arguments
#       ``epsilon``, ``sigma``, ``rc`` (cutoff), and ``smooth``.
#       Setting ``smooth=True`` applies a smooth cutoff to avoid
#       energy discontinuities.

# LJ parameters for argon (close to Rahman 1964 / White 1999 values)
eps_kB = 120.0       # K
sigma_Ar = 3.4       # Angstrom

# Convert epsilon from K to eV:  eps(eV) = eps(K) * kB(eV/K)
kB_eV = 8.617333e-5  # eV/K
epsilon_eV = eps_kB * kB_eV

# TODO: create a LennardJones calculator and assign it to atoms.calc
pass

# %%
# Let us verify the calculator works by computing the initial energy.

E0 = atoms.get_potential_energy()
print(f"Initial potential energy: {E0:.4f} eV")

# %%
# Initializing velocities
# ------------------------
#
# Before starting the MD we assign random velocities drawn from a
# Maxwell-Boltzmann distribution at the target temperature (85 K).

MaxwellBoltzmannDistribution(atoms, temperature_K=85)

# %%
# NVE vs NVT — the role of the thermostat
# -----------------------------------------
#
# Without a thermostat the simulation conserves total energy (NVE
# ensemble).  To maintain a constant *temperature* (NVT ensemble) we
# couple the system to a heat bath.
#
# **Langevin dynamics** adds two extra terms to the equation of motion:
#
# .. math::
#
#    m \, \mathbf{a} = \mathbf{F} - \gamma\,m\,\mathbf{v} + \boldsymbol{\xi}(t)
#
# The friction term :math:`-\gamma m \mathbf{v}` removes kinetic energy,
# while the random force :math:`\boldsymbol{\xi}(t)` (white noise) injects
# energy.  The balance between the two is set by the fluctuation-dissipation
# theorem so that, on average, the system reaches the desired temperature.
#
# Setting up the dynamics
# ------------------------
#
# .. seealso::
#
#    - `ASE Langevin dynamics <https://wiki.fysik.dtu.dk/ase/ase/md.html#ase.md.langevin.Langevin>`_
#    - `ASE MD module overview <https://wiki.fysik.dtu.dk/ase/ase/md.html>`_
#
# .. admonition:: Task 2
#
#    Create a ``Langevin`` dynamics object with a timestep of 5 fs,
#    a target temperature of 85 K, and a friction coefficient of
#    0.01 fs\ :sup:`-1`.
#
#    .. hint::
#
#       Use ``units.fs`` from ``ase.units`` to convert femtoseconds to
#       ASE internal time units.  The ``Langevin`` constructor takes
#       ``timestep``, ``temperature_K``, and ``friction``.

# TODO: create the Langevin dynamics object
dyn = None  # replace this

# %%
# Attaching output writers
# -------------------------
#
# We save the trajectory (positions, velocities, forces) and record the
# temperature and potential energy every 10 steps so we can analyse the
# simulation afterwards.
#
# .. admonition:: Task 3
#
#    Complete the ``record`` function to append the instantaneous temperature
#    and potential energy at each recording step.  Use ``atoms.get_temperature()``
#    and ``atoms.get_potential_energy()``.

traj = Trajectory("LJ_Argon_MD.traj", "w", atoms)
dyn.attach(traj.write, interval=10)

temperatures = []
energies = []


def record():
    pass  # TODO: append temperature and potential energy to the lists


dyn.attach(record, interval=10)

# %%
# Running the simulation
# -----------------------
#
# We run 1000 MD steps.  At each step the LJ forces are computed and the
# Velocity Verlet integrator (with Langevin thermostat) advances the
# positions and velocities.

n_steps = 1000
dyn.run(n_steps)
print(f"Ran {n_steps} steps.")

# %%
# Analysing the results
# ----------------------
#
# Let us plot the instantaneous temperature and potential energy as a
# function of simulation time.

time_ps = np.arange(len(temperatures)) * 10 * 5e-3  # steps * interval * dt(ps)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 5))

ax1.plot(time_ps, temperatures, "r-", linewidth=0.5)
ax1.axhline(85, color="k", linestyle="--", linewidth=0.8, label="Target: 85 K")
ax1.set_ylabel("Temperature (K)")
ax1.legend()

ax2.plot(time_ps, energies, "b-", linewidth=0.5)
ax2.set_ylabel("Potential energy (eV)")
ax2.set_xlabel("Time (ps)")

fig.suptitle("Liquid argon MD with Lennard-Jones potential")
fig.tight_layout()

# %%
# You can read the trajectory back for further analysis:
#
# .. code-block:: python
#
#     trajectory = ase.io.read("LJ_Argon_MD.traj", index=":")
#     # e.g. compute radial distribution function, diffusion coefficient, ...
#
