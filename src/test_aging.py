"""Eksperimen pelacakan umur & lead time (coupling komplementer Euler-Lagrange).

Menguji:
  1. Distribusi lead time produk yang selesai (histogram)
  2. Profil umur WIP terhadap tahap penyelesaian x
  3. Validasi hukum Little: E[WIP] = throughput_rate * E[lead time]

Hukum Little adalah uji fisis yang kuat: kalau tracking Lagrangian benar,
lead time rata-rata dikali laju throughput harus setara rata-rata WIP di
sistem pada kondisi tunak.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from adr_model import clearing_function
from aging import solve_aging_coupled, littles_law_check

OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def main():
    M = 200
    dx = 1.0 / M
    x = (np.arange(M) + 0.5) * dx
    mu = 1.0
    inflow = 0.5
    t_final = 8.0  # cukup panjang agar mencapai kondisi mendekati tunak

    # mulai dari sistem hampir kosong -> isi dari hulu
    rho0 = np.full(M, 0.05)

    (grid_x, times, rho_hist, lead_times, ages_final,
     ages_by_x, positions, birth_times) = solve_aging_coupled(
        rho0, t_final, dx, mu=mu, cfl=0.5, inflow=inflow, seed=1
    )

    rho_final = rho_hist[-1]
    wip_total = float(np.sum(rho_final) * dx)

    # laju throughput tunak = F(rho) di outflow (mendekati F(inflow) pada steady state)
    throughput_rate = clearing_function(np.array([rho_final[-1]]), mu)[0]

    N_pred, mean_LT = littles_law_check(lead_times, throughput_rate)

    print(f"[Aging] jumlah produk selesai (lead time tercatat) = {lead_times.size}")
    if lead_times.size:
        print(f"[Aging] lead time: mean={mean_LT:.4f}  std={lead_times.std():.4f}  "
              f"min={lead_times.min():.4f}  max={lead_times.max():.4f}")
    print(f"[Aging] WIP total akhir (E[WIP])         = {wip_total:.4f}")
    print(f"[Aging] throughput rate (F di outflow)   = {throughput_rate:.4f}")
    print(f"[Little] N_prediksi = throughput * LT    = {N_pred:.4f}")
    print(f"[Little] rasio N_pred / WIP_aktual       = {N_pred / wip_total:.3f}  (target ~1.0)")

    # analitis: pada rho seragam = inflow, kecepatan v = mu/(1+inflow),
    # lead time = 1 / v = (1+inflow)/mu
    LT_analitis = (1.0 + inflow) / mu
    print(f"[Aging] lead time analitis (rho~inflow)  = {LT_analitis:.4f}")

    # ---- plot ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    # (a) profil kerapatan akhir
    axes[0].plot(x, rho_final, color="steelblue")
    axes[0].axhline(inflow, color="gray", ls="--", alpha=0.6, label=f"inflow={inflow}")
    axes[0].set_xlabel("x (derajat penyelesaian)")
    axes[0].set_ylabel(r"$\rho$ (WIP)")
    axes[0].set_title(f"Profil WIP pada t={t_final:.0f} (mendekati tunak)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    # (b) distribusi lead time
    if lead_times.size:
        axes[1].hist(lead_times, bins=40, color="darkorange", edgecolor="white")
        axes[1].axvline(mean_LT, color="red", ls="--", label=f"mean={mean_LT:.3f}")
        axes[1].axvline(LT_analitis, color="black", ls=":", label=f"analitis={LT_analitis:.3f}")
        axes[1].legend(fontsize=8)
    axes[1].set_xlabel("lead time (satuan waktu)")
    axes[1].set_ylabel("jumlah produk")
    axes[1].set_title("Distribusi lead time produk selesai")
    axes[1].grid(alpha=0.3)

    # (c) profil umur WIP vs x
    valid = ~np.isnan(ages_by_x)
    axes[2].plot(x[valid], ages_by_x[valid], color="green")
    axes[2].set_xlabel("x (derajat penyelesaian)")
    axes[2].set_ylabel("umur rata-rata WIP (waktu)")
    axes[2].set_title("Profil umur WIP terhadap tahap")
    axes[2].grid(alpha=0.3)

    fig.suptitle("Pelacakan umur & lead time (coupling komplementer Euler-Lagrange)",
                 fontsize=12)
    fig.tight_layout()
    path = os.path.join(OUT, "06_aging_leadtime.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[Aging] Plot disimpan: {path}")


if __name__ == "__main__":
    main()
