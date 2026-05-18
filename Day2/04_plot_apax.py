#!/usr/bin/python

""".. _plot_apax:

Exercise 2.3: Evaluating the trained Apax model
=================================================

In this tutorial we evaluate how well the trained Apax model reproduces
the DFT reference data by creating **energy parity plots**. In a parity
plot we scatter the DFT reference energy (x-axis) against the
Apax-predicted energy (y-axis) for every frame. If the model were
perfect, all points would fall exactly on the diagonal line.

Let us start by importing our Python modules.
"""

import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from apax.md import ASECalculator
from pathlib import Path

# %%
# Set the project path and load the trained model.

PROJECT_PATH = Path("../../../solutions")
model_dir = PROJECT_PATH / "apax" / "Argon_Simulation"
calc = ASECalculator(model_dir=model_dir)

# %%
# Helper functions
# ----------------
#
# We define a small helper to compute the root mean square error (RMSE)
# between reference and predicted values.


def rms_dict(x_ref, x_pred):
    """Compute RMSE and standard deviation of errors."""
    x_ref = np.array(x_ref)
    x_pred = np.array(x_pred)
    error = x_ref - x_pred
    rmse = np.sqrt(np.mean(error**2))
    std = np.std(error)
    return {"rmse": rmse, "std": std}


# %%
# Next, we write a function that creates an energy parity plot for a
# given data file.
#
# .. admonition:: Task 1
#
#    Implement the ``energy_parity_plot`` function below. For each frame
#    in the data file you need to:
#
#    1. Read the DFT reference energy from ``frame.info["energy"]``
#    2. Attach the Apax calculator and compute the predicted energy
#    3. Normalize both energies per atom
#    4. Create a scatter plot and add the diagonal reference line


def energy_parity_plot(data_file, ax, title="Energy parity plot"):
    """Create an energy parity plot comparing DFT and Apax energies."""
    pass  # TODO: implement the energy parity plot


# %%
# Finally, we create the parity plots for the training and test data
# side by side.

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

energy_parity_plot(
    PROJECT_PATH / "apax/train.xyz",
    axes[0],
    title="Energy on training data",
)
energy_parity_plot(
    PROJECT_PATH / "apax/test.xyz",
    axes[1],
    title="Energy on test data",
)

fig.tight_layout()
fig.savefig("apax_energy_parity.png", dpi=150)
plt.show()
