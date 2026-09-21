"""Eksperimen: efek bottleneck kapasitas terhadap WIP, antrian, dan lead time.

Bandingkan dua skenario:
  A. Tanpa bottleneck (kapasitas seragam mu=1.0)
  B. Dengan bottleneck (kapasitas turun ke 0.35 pada x in [0.55, 0.65])

Diharapkan pada skenario B:
  - WIP menumpuk di HULU bottleneck (antrian merambat mundur -- back-pressure)
  - lead time membengkak & distribusinya melebar
  - profil umur naik lebih curam di depan bottleneck
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from bottleneck import solve_bottleneck_aging

OUT = os.path.join(os.path.dirname(__file__), "..", "results")


def summarize(tag, res):
    lt = res["lead_times"]
    wip = float(np.sum(res["rho_final"]) * res["dx"])
    if lt.size:
        print(f"[{tag}] produk selesai={lt.size}  lead time: mean={lt.mean():.3f} "
              f"std={lt.std():.3f} p95={np.percentile(lt,95):.3f}")
    else:
        print(f"[{tag}] belum ada produk selesai")
    print(f"[{tag}] WIP total akhir = {wip:.3f}")
    return lt, wip


def main():
    common = dict(M=200, t_final=25.0, q_in=0.5, rho_max=8.0, cfl=0.4, seed=1)

    # A. tanpa bottleneck
    res_a = solve_bottleneck_aging(mu_normal=1.0, mu_bottleneck=1.0, **common)
    # B. dengan bottleneck
    res_b = solve_bottleneck_aging(mu_normal=1.0, mu_bottleneck=0.35,
                                   x_lo=0.55, x_hi=0.65, **common)

    print("=" * 60)
    lt_a, wip_a = summarize("Tanpa bottleneck", res_a)
    lt_b, wip_b = summarize("Dengan bottleneck", res_b)
    print("=" * 60)
    if lt_a.size and lt_b.size:
        print(f"[Efek] lead time mean: {lt_a.mean():.3f} -> {lt_b.mean():.3f} "
              f"(naik {lt_b.mean()/lt_a.mean():.2f}x)")
    print(f"[Efek] WIP total: {wip_a:.3f} -> {wip_b:.3f} (naik {wip_b/wip_a:.2f}x)")

    x = res_a["grid_x"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # (a) profil WIP akhir
    ax = axes[0, 0]
    ax.plot(x, res_a["rho_final"], color="steelblue", label="tanpa bottleneck")
    ax.plot(x, res_b["rho_final"], color="crimson", label="dengan bottleneck")
    ax.axvspan(res_b["x_lo"], res_b["x_hi"], color="gray", alpha=0.2,
               label="zona bottleneck")
    ax.set_xlabel("x (derajat penyelesaian)")
    ax.set_ylabel(r"$\rho$ (WIP)")
    ax.set_title("Profil WIP akhir: buffer hulu jenuh akibat back-pressure")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (b) profil kapasitas
    ax = axes[0, 1]
    ax.plot(x, res_b["mu"], color="black")
    ax.axvspan(res_b["x_lo"], res_b["x_hi"], color="gray", alpha=0.2)
    ax.set_xlabel("x (derajat penyelesaian)")
    ax.set_ylabel(r"$\mu(x)$ kapasitas")
    ax.set_title("Profil kapasitas per tahap")
    ax.grid(alpha=0.3)

    # (c) distribusi lead time
    ax = axes[1, 0]
    if lt_a.size:
        ax.hist(lt_a, bins=40, alpha=0.6, color="steelblue", label="tanpa bottleneck")
    if lt_b.size:
        ax.hist(lt_b, bins=40, alpha=0.6, color="crimson", label="dengan bottleneck")
    ax.set_yscale("log")
    ax.set_xlabel("lead time (satuan waktu)")
    ax.set_ylabel("jumlah produk (skala log)")
    ax.set_title("Distribusi lead time: membengkak & melebar")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # (d) profil umur WIP vs x
    ax = axes[1, 1]
    for res, c, lab in [(res_a, "steelblue", "tanpa bottleneck"),
                        (res_b, "crimson", "dengan bottleneck")]:
        a = res["ages_by_x"]
        v = ~np.isnan(a)
        ax.plot(x[v], a[v], color=c, label=lab)
    ax.axvspan(res_b["x_lo"], res_b["x_hi"], color="gray", alpha=0.2)
    ax.set_xlabel("x (derajat penyelesaian)")
    ax.set_ylabel("umur rata-rata WIP (waktu)")
    ax.set_title("Profil umur WIP: melonjak di depan bottleneck")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle("Dampak bottleneck kapasitas terhadap WIP, antrian & lead time",
                 fontsize=13)
    fig.tight_layout()
    path = os.path.join(OUT, "07_bottleneck_effect.png")
    fig.savefig(path, dpi=125)
    plt.close(fig)
    print(f"[Bottleneck] Plot disimpan: {path}")


if __name__ == "__main__":
    main()
