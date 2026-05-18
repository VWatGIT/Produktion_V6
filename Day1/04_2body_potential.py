#!/usr/bin/python

""".. _2_body_potential:

Exercise 1.3: Calculating the 2 body potential of Argon using ASE and CP2K
===========================================================================

A **2-body potential** (also called a *pair potential*) describes the interaction
energy between two atoms as a function of their separation :math:`r` only:

.. math::

   V_\\text{2-body}(r) = V(r_{12})

The total potential energy of a many-body system is then approximated as a sum
over all distinct pairs:

.. math::

   E = \\sum_{i < j} V(r_{ij})

This is the class of potential that the Lennard-Jones model belongs to — each
atom pair contributes independently, and three-body (and higher) correlations
are neglected.

In this tutorial we compute the 2-body potential of argon *from first principles*
using DFT, and compare it with the Lennard-Jones curve fitted in Exercise 1.1.
Because we are sampling only two atoms, the calculation is cheap enough to be
done for many separations.

We use ASE (Atomic Simulation Environment) to set up the system.
CP2K is called by ASE to compute the DFT energy at each geometry.

First, we import our Python modules.
"""

# uncomment the following line if running a jupyter notebook
# %matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms, units
from ase.calculators.cp2k import CP2K

# %%
# We also import pint, which will help us to keep track of units.

import pint

ureg = pint.UnitRegistry()

# %%
# Next, we create an Atoms object with two argon atoms placed a distance apart.

distance = 3.3 * ureg.angstrom
two_argon_atoms = Atoms("Ar2", [[0, 0, 0], [0, 0, distance.magnitude]])

# %%
# We set the cubic simulation box size to something much larger than Ar-Ar distance
# and apply periodic boundary conditions in all directions.

two_argon_atoms.center(vacuum=3)
two_argon_atoms.pbc = [1, 1, 1]

# %%
# We will use CP2K as a calculator for our Atoms object.
# First, we specify calculation settings, you will not have to worry about these

inp = """
# Parameters for force calculation.
&FORCE_EVAL
    # Define the DFT parameters

    &DFT
        &SCF
            &OT
                MINIMIZER DIIS
                PRECONDITIONER FULL_SINGLE_INVERSE
            &END OT
            &OUTER_SCF
                MAX_SCF 100
                EPS_SCF 1.0E-6
            &END OUTER_SCF
        &END SCF

    &END DFT
&END FORCE_EVAL
"""

# %%
# Next, we specify our basis set, pseudo potential and exchange correlation functional.
#
# .. seealso::
#
#    - `ASE CP2K calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/cp2k.html>`_
#    - `CP2K BASIS_MOLOPT basis sets <https://github.com/cp2k/cp2k/blob/master/data/BASIS_MOLOPT>`_
#
# .. admonition:: Task 1
#
#    Create the CP2K calculator below and assign it to ``two_argon_atoms.calc``.
#    Use ``basis_set="DZVP-MOLOPT-SR-GTH"``, ``pseudo_potential="GTH-PBE-q8"``,
#    ``potential_file="GTH_POTENTIALS"``, and ``xc="PBE"``.
#
#    The ``cutoff`` parameter controls the plane-wave grid resolution (in Ry).
#    A higher cutoff gives more accurate results but is more expensive.
#    Try different values (e.g. 200, 400, 600 Ry) and observe how the total
#    energy changes.  At what cutoff is the energy converged to within ~1 meV?
#
#    .. hint::
#
#       Pass the cutoff via the ``cutoff`` keyword of the ``CP2K`` constructor.
#       Note that CP2K uses Rydberg (Ry) as the cutoff unit.

# TODO: create the CP2K calculator and assign it to two_argon_atoms.calc
# Remember to include a cutoff parameter (try e.g. 400 Ry)
two_argon_atoms.calc = CP2K(inp=inp,
                            basis_set="DZVP-MOLOPT-SR-GTH",
                            pseudo_potential="GTH-PBE-q8",
                            potential_file="GTH_POTENTIALS",
                            xc="PBE",
                            cutoff=400 * units.Rydberg,
                            command="cp2k.psmp -s",
                            )


# %%
# We can now run our first DFT calculation with

E = two_argon_atoms.get_potential_energy()
print(E)

# %%
# This will return us a single potential energy for the specified interatomic distance.
# Next, we want to sample a region between 3.3 and 6.0 Angstrom and get the potential energy for each distance

distances = np.linspace(3.3, 6.0, 20)
energies = np.zeros(distances.shape)

for i in range(len(distances)):
    two_argon_atoms.set_positions([[0, 0, 0], [0, 0, distances[i]]])
    two_argon_atoms.center(vacuum=3)
    print(distances[i])
    print(two_argon_atoms.get_positions())

    E = two_argon_atoms.get_potential_energy()
    print(E)
    energies[i] = E

# %%
# Finally, we plot the energy as a function of the interatomic distance.
#
# .. admonition:: Task 2
#
#    Plot the DFT 2-body potential below.  Then overlay the Lennard-Jones
#    potential you fitted in Exercise 1.1 (with :math:`\varepsilon/k_B = 120\,\text{K}`
#    and :math:`\sigma = 3.4\,\text{\AA}`) on the same axes so you can
#    compare the two.
#
#    Hint: use the ``lj_potential`` formula
#    :math:`V(r) = 4\varepsilon[(\sigma/r)^{12} - (\sigma/r)^{6}]`
#    and convert from Joules to eV.

# TODO: plot the DFT energies and overlay the LJ potential from Day 1
