#!/usr/bin/python

""".. _rdf:

Exercise 3.2: Calculating the radial distribution function
************************************************************

In this tutorial we learn how to analyse molecular dynamics trajectories.
We use MDAnalysis to read our ASE trajectory from the previous exercise
and calculate the radial distribution function (RDF).

`Wikipedia page on RDF <https://en.wikipedia.org/wiki/Radial_distribution_function>`_

`MDAnalysis documentation on RDF <https://docs.mdanalysis.org/stable/documentation_pages/analysis/rdf.html>`_

We first import our Python modules.
"""

import matplotlib.pyplot as plt
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

# %%
# Then we define our project path. Replace the path with your own project path.

PROJECT_PATH = Path("../../../solutions/")

# %%
# Loading the trajectory
# ----------------------
#
# We load the ASE trajectory produced by the Apax MD simulation.
# MDAnalysis cannot read ``.traj`` files directly, so we first convert
# to extxyz using ASE, then load with MDAnalysis.
#
# .. admonition:: Task 1
#
#    Load the trajectory ``Argon_ML_85K.traj`` with ``ase.io.read``,
#    write it to a temporary extxyz file, then create an MDAnalysis
#    ``Universe`` from it. Set the box dimensions to
#    17.0742 Å in each direction.

universe = None  # TODO: load trajectory, convert to xyz, create Universe

# %%
# The Universe object contains the atomic positions for each timestep.
# Note that the xyz file does not contain box dimension information that
# MDAnalysis can read, so we must set it ourselves.

box_l = 17.0742
universe.dimensions = [box_l, box_l, box_l, 90, 90, 90]

# %%
# Let us check how many frames we have loaded.

print(f"loaded {len(universe.trajectory)} frames")

# %%
# Computing the RDF
# ------------------
#
# We now run a radial distribution function analysis using InterRDF.

rdf = InterRDF(
    universe.atoms, universe.atoms,
    n_bins=100,
    range=(1.0, box_l / 2),
)

# %%
# We run the analysis.

rdf.run()

# %%
# Plotting the results
# ---------------------

plt.plot(rdf.results.bins, rdf.results.rdf)
plt.xlabel(r"$r$ ($\AA$)")
plt.ylabel("g(r)")
plt.title("Radial distribution function — Apax MD at 85 K")
plt.tight_layout()

# %%
# Save the figure.

plt.savefig("rdf_apax_85K.png", dpi=300)
