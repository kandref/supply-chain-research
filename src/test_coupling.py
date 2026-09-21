"""Uji coupling multi-skala Euler-Lagrange (splicing) pada bottleneck (shock).

Verifikasi:
  - massa hasil coupling lestari (dekat Eulerian & Lagrangian)
  - posisi shock cocok dengan Rankine-Hugoniot
  - coupling TIDAK merusak solusi; profil masuk akal di seluruh domain
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adr_model import clearing_function
from coupling import solve_coupled

OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def main():
    M = 400
    dx = 1.0 / M
    x = (np.arange(M) + 0.5) * dx
    t_final = 0.6
    mu = 1.0

    rho_l, rho_r = 0.3, 2.5
    rho0 = np.where(x < 0.5, rho_l, rho_r)

    grid_x, rho_eul, rho_lag, rho_cpl, window = solve_coupled(
        rho0, t_final, dx, mu=mu,
        cfl_eul=0.5, cfl_lag=0.4,
        window_half_width=30, blend_width=8,
        parcels_per_unit_mass=15000, bandwidth_factor=1.5,
    )

    fl = clearing_function(np.array([rho_l]), mu)[0]
    fr = clearing_function(np.array([rho_r]), mu)[0]
    s_rh = (fr - fl) / (rho_r - rho_l)
    shock_x_ana = 0.5 + s_rh * t_final
    rho_mid = 0.5 * (rho_l + rho_r)

    def shock_pos(r):
        return x[np.argmin(np.abs(r - rho_mid))]

    mass = lambda r: float(np.sum(r) * dx)
    err_l1 = np.mean(np.abs(rho_eul - rho_cpl))

    print(f"[Coupled] Rankine-Hugoniot analitis : x_shock = {shock_x_ana:.4f}")
    print(f"[Coupled] shock Eulerian            : x = {shock_pos(rho_eul):.4f}")
    print(f"[Coupled] shock Lagrangian          : x = {shock_pos(rho_lag):.4f}")
    print(f"[Coupled] shock Coupled             : x = {shock_pos(rho_cpl):.4f}")
    print(f"[Coupled] massa: Eul={mass(rho_eul):.5f}  Lag={mass(rho_lag):.5f}  Cpl={mass(rho_cpl):.5f}")
    print(f"[Coupled] galat ||Eul - Coupled||_L1 = {err_l1:.4e}")
    print(f"[Coupled] jendela Lagrangian (i_lo, i_hi) = {window} -> x in [{x[window[0]]:.3f}, {x[window[1]]:.3f}]")

    # ukur ketajaman shock: lebar transisi 10%-90%
    def shock_width(r):
        lo, hi = rho_l + 0.1 * (rho_r - rho_l), rho_l + 0.9 * (rho_r - rho_l)
        i1 = np.argmin(np.abs(r - lo))
        i2 = np.argmin(np.abs(r - hi))
        return abs(x[i2] - x[i1])

    print(f"[Coupled] lebar transisi shock: Eul={shock_width(rho_eul):.4f}  Cpl={shock_width(rho_cpl):.4f}")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, rho0, ":", color="gray", alpha=0.5, label=r"$\rho(x,0)$")
    ax.plot(x, rho_eul, "-", color="steelblue", lw=2, label="Eulerian murni")
    ax.plot(x, rho_lag, "-.", color="green", lw=1, alpha=0.7, label="Lagrangian penuh")
    ax.plot(x, rho_cpl, "--", color="crimson", lw=1.8, label="Coupled (splicing)")
    ax.axvspan(x[window[0]], x[window[1]], color="orange", alpha=0.12,
               label="jendela Lagrangian")
    ax.axvline(shock_x_ana, color="black", ls=":", alpha=0.7,
               label=f"shock RH = {shock_x_ana:.3f}")
    ax.set_xlabel("x (derajat penyelesaian)")
    ax.set_ylabel(r"$\rho$ (WIP)")
    ax.set_title(f"Coupling multi-skala pada bottleneck, t = {t_final}")
    ax.legend(fontsize=8.5, loc="center left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    path = os.path.join(OUT, "05_coupling.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[Coupled] Plot disimpan: {path}")


if __name__ == "__main__":
    main()
