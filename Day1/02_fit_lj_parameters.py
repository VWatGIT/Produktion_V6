#!/usr/bin/python

""".. _fit_lj_parameters:

Exercise 1.1: Fitting Lennard-Jones parameters for argon
=========================================================

In this tutorial we learn how to determine Lennard-Jones (LJ) potential
parameters for argon by fitting to the second virial coefficient *B(T)*.
Along the way we derive the statistical-mechanical expression for *B(T)*
and compute the integral numerically.

First, we import our Python modules.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import least_squares

# %%
# The Lennard-Jones potential
# ---------------------------
#
# The 12-6 Lennard-Jones potential describes the interaction between two
# neutral atoms as a function of their separation *r*:
#
# .. math::
#
#    V(r) = 4\varepsilon \left[\left(\frac{\sigma}{r}\right)^{12}
#    - \left(\frac{\sigma}{r}\right)^{6}\right]
#
# Here :math:`\varepsilon` is the depth of the potential well and
# :math:`\sigma` is the distance at which the potential crosses zero.
# The :math:`r^{-6}` term models the attractive van der Waals (London
# dispersion) interaction, while the :math:`r^{-12}` term is a steep
# repulsive wall at short range.
#
# .. admonition:: Task 1
#
#    Implement the ``lj_potential`` function below. Given the interatomic
#    distance ``r``, well depth ``epsilon``, and zero-crossing distance
#    ``sigma``, return the LJ potential energy :math:`V(r)`.


def lj_potential(r, epsilon, sigma):
    """Lennard-Jones pair potential.

    Parameters
    ----------
    r : array_like
        Interatomic distance (m).
    epsilon : float
        Well depth (J).
    sigma : float
        Zero-crossing distance (m).

    Returns
    -------
    V : array_like
        Potential energy (J).
    """
    pass  # TODO: implement the LJ potential formula


# %%
# Let us plot the LJ potential for some typical argon parameters.

K_B = 1.380649e-23     # Boltzmann constant, J/K
N_A = 6.02214076e23    # Avogadro's number

# Starting guess — close to well-known literature values
sigma_guess = 3.4e-10   # m
eps_guess = 120.0 * K_B  # J

r = np.linspace(2.8e-10, 8.0e-10, 500)
V = lj_potential(r, eps_guess, sigma_guess)

plt.figure()
plt.plot(r * 1e10, V / K_B, "k-")
plt.axhline(0, color="grey", linewidth=0.5)
plt.xlabel(r"$r$ ($\AA$)")
plt.ylabel(r"$V(r)/k_B$ (K)")
plt.title("Lennard-Jones potential for argon")
plt.ylim(-200, 400)
plt.tight_layout()

# %%
# The second virial coefficient
# ---------------------------------------
#
# Deriving B₂(T) from statistical mechanics
# ------------------------------------------
#
# The equation of state for a real gas can be expanded in powers of the
# number density :math:`\rho = N/V`:
#
# .. math::
#
#    \frac{P}{k_B T} = \rho + B_2(T)\,\rho^2 + B_3(T)\,\rho^3 + \cdots
#
# The second virial coefficient :math:`B_2(T)` captures the leading-order
# correction to ideal-gas behaviour and depends only on the pair potential.
#
# For a spherically symmetric pair potential :math:`V(r)`, statistical
# mechanics gives:
#
# .. math::
#
#    B_2(T) = -2\pi N_A \int_0^{\infty}
#    \left[\exp\!\left(-\frac{V(r)}{k_B T}\right) - 1\right] r^2 \, dr
#
# The factor :math:`N_A` converts from per-molecule to per-mole (so that
# :math:`B_2` has units of cm³/mol, matching the experimental data).
#
# For the LJ potential this integral cannot be evaluated analytically, so
# we compute it numerically.
#
# .. admonition:: Task 2
#
#    Implement the ``second_virial`` function below. Use reduced units
#    :math:`r^* = r/\sigma`, :math:`T^* = k_B T/\varepsilon` for numerical
#    stability. The integral becomes:
#
#    .. math::
#
#       B_2 = -2\pi N_A \sigma^3 \int_0^{\infty}
#       \left[\exp\!\left(-\frac{V^*(r^*)}{T^*}\right) - 1\right] r^{*2} \, dr^*
#
#    where :math:`V^*(r^*) = 4(r^{*-12} - r^{*-6})`.
#    Use ``scipy.integrate.quad`` for the numerical integration.
#    Return the result in cm³/mol.


def second_virial(T, epsilon, sigma):
    """Compute B2(T) for the LJ potential by numerical integration.

    Parameters
    ----------
    T : float
        Temperature (K).
    epsilon, sigma : float
        LJ parameters (J, m).

    Returns
    -------
    B2 : float
        Second virial coefficient (cm³/mol).
    """
    pass  # TODO: implement the second virial coefficient calculation


# %%
# Experimental data
# -----------------
#
# We use the smoothed second virial coefficients from Dymond & Smith,
# *The Virial Coefficients of Pure Gases and Mixtures* (Oxford, 1980).
# Units are cm³/mol.

T_B_exp = np.array([81, 85, 90, 95, 100, 110, 125, 150,
                     200, 250, 300, 400, 500, 600, 700, 800, 900, 1000],
                    dtype=float)
B_exp = np.array([-276, -251, -225, -202.5, -183.5, -154.5, -123.0, -86.2,
                   -47.4, -27.9, -15.5, -1.0, 7.0, 12.0, 15.0, 17.7, 20.0, 22.0])

# %%
# Fitting ε and σ to B(T)
# ------------------------
#
# We minimize the sum of squared residuals between the computed and
# experimental :math:`B_2(T)` values using ``scipy.optimize.least_squares``.
#
# .. admonition:: Task 3
#
#    Implement the ``residuals_B`` function and run the least-squares fit.
#    The function should take ``params = [sigma, eps_over_kB]`` and return
#    the vector of residuals ``B_calc - B_exp``.


def residuals_B(params):
    sigma, eps_over_kB = params
    epsilon = eps_over_kB * K_B
    pass  # TODO: compute B_calc for each temperature and return residuals


result_B = least_squares(residuals_B,
                          x0=[3.4e-10, 120.0],
                          bounds=([2.5e-10, 50.0], [5.0e-10, 250.0]))
sigma_B, eps_kB_B = result_B.x
print(f"Fit from B(T):  sigma = {sigma_B*1e10:.4f} A,  epsilon/kB = {eps_kB_B:.2f} K")

# %%
# Let us compare the fit with the experimental data.

T_fine = np.linspace(80, 1000, 200)
B_fit = np.array([second_virial(T, eps_kB_B * K_B, sigma_B) for T in T_fine])

plt.figure()
plt.plot(T_B_exp, B_exp, "ko", label="Experiment (Dymond & Smith)")
plt.plot(T_fine, B_fit, "r-", label="LJ fit")
plt.axhline(0, color="grey", linewidth=0.5)
plt.xlabel("$T$ (K)")
plt.ylabel("$B_2$ (cm$^3$/mol)")
plt.legend()
plt.title("Second virial coefficient of argon")
plt.tight_layout()

# %%
# Comparing with literature values
# ---------------------------------
#
# Let us compare the LJ parameters obtained from the fit with
# literature values.

print("\n--- Comparison of LJ parameters for argon ---")
print(f"{'Source':<30s} {'eps/kB (K)':>12s} {'sigma (A)':>12s}")
print("-" * 56)
print(f"{'Fit from B(T)':<30s} {eps_kB_B:12.2f} {sigma_B*1e10:12.4f}")
print(f"{'Rahman (1964)':<30s} {'120.0':>12s} {'3.4000':>12s}")
print(f"{'White (1999)':<30s} {'125.7':>12s} {'3.3450':>12s}")

# %%
# The fact that a simple two-parameter potential can reproduce the
# experimental virial coefficients to within a few percent is a testament
# to the usefulness of the LJ model — but also highlights its limitations.
#
