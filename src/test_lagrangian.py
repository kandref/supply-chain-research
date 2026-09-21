"""Uji Lagrangian vs Eulerian pada kondisi awal mulus.

Kalau dua deskripsi konsisten, keduanya harus menghasilkan profil akhir yang
serupa (dalam batas galat rekonstruksi KDE + galat Godunov orde-1).
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adr_model import solve_adr
from lagrangian import solve_lagrangian

OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def main():
    M = 200
    dx = 1.0 / M
    x = (np.arange(M) + 0.5) * dx
    t_final = 0.3
    mu = 1.0

    # kondisi awal mulus, positif, tidak terlalu tinggi
    rho0 = 1.0 + 0.5 * np.sin(2 * np.pi * x)

    # solver Eulerian (referensi)
    _, hist_eul = solve_adr(rho0, t_final, dx, mu=mu, cfl=0.5, save_every=10_000)
    rho_eul = hist_eul[-1]

    # solver Lagrangian
    _, rho_lag = solve_lagrangian(
        rho0, t_final, dx, mu=mu, cfl=0.4, parcels_per_unit_mass=8000, bandwidth_factor=2.0
    )

    err_l1 = np.mean(np.abs(rho_eul - rho_lag))
    err_linf = np.max(np.abs(rho_eul - rho_lag))
    mass_eul = float(np.sum(rho_eul) * dx)
    mass_lag = float(np.sum(rho_lag) * dx)
    mass_init = float(np.sum(rho0) * dx)

    print(f"[X-check] galat ||Eul - Lag||_L1  = {err_l1:.4e}")
    print(f"[X-check] galat ||Eul - Lag||_Linf = {err_linf:.4e}")
    print(f"[X-check] massa awal      = {mass_init:.5f}")
    print(f"[X-check] massa Eulerian  = {mass_eul:.5f}")
    print(f"[X-check] massa Lagrangian= {mass_lag:.5f}")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, rho0, "--", color="gray", label="rho(x, 0)")
    ax.plot(x, rho_eul, "-", color="steelblue", label="Eulerian (Godunov)")
    ax.plot(x, rho_lag, "--", color="darkorange", label="Lagrangian (partikel + KDE)")
    ax.set_xlabel("x (derajat penyelesaian)")
    ax.set_ylabel(r"$\rho$ (WIP)")
    ax.set_title(f"Cross-check Eulerian vs Lagrangian, t = {t_final}")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = os.path.join(OUT, "04_xcheck_eul_lag.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[X-check] Plot disimpan: {path}")


if __name__ == "__main__":
    main()
