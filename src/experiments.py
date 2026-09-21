"""Eksperimen numerik verifikasi solver ADR.

Menjalankan tiga kasus utama:
  1. Riemann problem: pembentukan bottleneck (shock) & pelepasan (rarefaction)
  2. Uji konvergensi terhadap ukuran grid
  3. Kuantifikasi ketidakpastian Monte-Carlo pada throughput
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib

matplotlib.use("Agg")  # backend headless
import matplotlib.pyplot as plt

from adr_model import (
    clearing_function,
    solve_adr,
    step_euler,
    cfl_dt,
    total_mass,
)

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT, exist_ok=True)


# --------------------------------------------------------------------------
# Eksperimen 1: Riemann problem -- bottleneck (shock) & pelepasan (rarefaction)
# --------------------------------------------------------------------------
def run_riemann(mu: float = 1.0, M: int = 400, t_final: float = 0.6):
    dx = 1.0 / M
    x = (np.arange(M) + 0.5) * dx

    # Kasus A: shock -- kerapatan rendah di hulu, tinggi di hilir
    # (analog: hilir tersumbat -> WIP menumpuk merambat mundur)
    rho0_shock = np.where(x < 0.5, 0.3, 2.5)
    t_s, h_s = solve_adr(rho0_shock, t_final, dx, mu=mu, cfl=0.8, save_every=10)

    # Kasus B: rarefaction -- kerapatan tinggi di hulu, rendah di hilir
    # (analog: bottleneck baru saja dilepaskan)
    rho0_rare = np.where(x < 0.5, 2.5, 0.3)
    t_r, h_r = solve_adr(rho0_rare, t_final, dx, mu=mu, cfl=0.8, save_every=10)

    # Plot: dua panel
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    panels = [
        ("Bottleneck (shock) merambat mundur", x, h_s, t_s),
        ("Pelepasan bottleneck (rarefaction)", x, h_r, t_r),
    ]
    for ax, (title, x_, hist, times) in zip(axes, panels):
        # tampilkan beberapa snapshot
        idxs = np.linspace(0, len(hist) - 1, 5).astype(int)
        for i in idxs:
            ax.plot(x_, hist[i], label=f"t = {times[i]:.2f}")
        ax.set_xlabel("x (derajat penyelesaian)")
        ax.set_ylabel(r"$\rho$ (WIP)")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.suptitle("Kasus kanonik: Riemann problem pada model ADR", fontsize=12)
    fig.tight_layout()
    path = os.path.join(OUT, "01_riemann.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)

    # Space-time diagram bottleneck
    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(
        h_s,
        aspect="auto",
        origin="lower",
        extent=[0, 1, 0, t_s[-1]],
        cmap="magma",
    )
    ax.set_xlabel("x (derajat penyelesaian)")
    ax.set_ylabel("t (waktu)")
    ax.set_title("Space-time diagram: bottleneck (shock) merambat mundur")
    fig.colorbar(im, ax=ax, label=r"$\rho$ (WIP)")
    fig.tight_layout()
    st_path = os.path.join(OUT, "01_riemann_spacetime.png")
    fig.savefig(st_path, dpi=130)
    plt.close(fig)

    # Rankine-Hugoniot check
    rho_l, rho_r = 0.3, 2.5
    fl = clearing_function(np.array([rho_l]), mu)[0]
    fr = clearing_function(np.array([rho_r]), mu)[0]
    s_rh = (fr - fl) / (rho_r - rho_l)

    # ukur kecepatan shock dari simulasi: cari titik tengah lonjakan pada t akhir
    rho_end = h_s[-1]
    rho_mid = 0.5 * (rho_l + rho_r)
    shock_x_num = x[np.argmin(np.abs(rho_end - rho_mid))]
    shock_x_ana = 0.5 + s_rh * t_s[-1]

    print(f"[Riemann] Rankine-Hugoniot s_analitis = {s_rh:+.4f}")
    print(f"[Riemann] Posisi shock analitis pada t={t_s[-1]:.3f}: x = {shock_x_ana:.4f}")
    print(f"[Riemann] Posisi shock numerik   pada t={t_s[-1]:.3f}: x = {shock_x_num:.4f}")
    print(f"[Riemann] Selisih |x_num - x_ana| = {abs(shock_x_num - shock_x_ana):.4f}")
    print(f"[Riemann] Plot disimpan: {path} & {st_path}")
    return s_rh, shock_x_num, shock_x_ana


# --------------------------------------------------------------------------
# Eksperimen 2: Uji konvergensi
# --------------------------------------------------------------------------
def run_convergence(mu: float = 1.0, t_final: float = 0.3):
    """Uji konvergensi terhadap solusi referensi (grid halus).

    Kondisi awal mulus (sinusoidal) agar solusi belum diskontinu -- ini yang
    membuat orde konvergensi teoretis (~1 untuk Godunov orde satu) terukur
    dengan bersih.
    """
    Ms = [50, 100, 200, 400]
    M_ref = 1600
    dx_ref = 1.0 / M_ref
    x_ref = (np.arange(M_ref) + 0.5) * dx_ref
    rho0_ref = 1.0 + 0.5 * np.sin(2 * np.pi * x_ref)
    _, hist_ref = solve_adr(rho0_ref, t_final, dx_ref, mu=mu, cfl=0.5, save_every=10_000)
    rho_ref = hist_ref[-1]

    errors_l1 = []
    for M in Ms:
        dx = 1.0 / M
        x = (np.arange(M) + 0.5) * dx
        rho0 = 1.0 + 0.5 * np.sin(2 * np.pi * x)
        _, hist = solve_adr(rho0, t_final, dx, mu=mu, cfl=0.5, save_every=10_000)
        rho_end = hist[-1]
        # bandingkan lewat rata-rata blok referensi
        block = M_ref // M
        rho_ref_avg = rho_ref.reshape(M, block).mean(axis=1)
        err = np.mean(np.abs(rho_end - rho_ref_avg))  # norma L1 diskret
        errors_l1.append(err)
        print(f"[Konvergensi] M={M:4d}  dx={dx:.5f}  ||err||_1 = {err:.5e}")

    # Estimasi orde konvergensi dari pasangan berturut-turut
    orders = []
    for i in range(1, len(Ms)):
        p = np.log(errors_l1[i - 1] / errors_l1[i]) / np.log(Ms[i] / Ms[i - 1])
        orders.append(p)
        print(f"[Konvergensi] M={Ms[i-1]}->M={Ms[i]}  orde empiris p ~ {p:.3f}")

    fig, ax = plt.subplots(figsize=(6, 4.5))
    dxs = [1.0 / M for M in Ms]
    ax.loglog(dxs, errors_l1, "o-", label="galat L1 numerik")
    # referensi orde 1
    ref = errors_l1[0] * (np.array(dxs) / dxs[0]) ** 1.0
    ax.loglog(dxs, ref, "--", color="gray", label="orde 1 (referensi)")
    ax.set_xlabel(r"$\Delta x$")
    ax.set_ylabel("galat L1")
    ax.set_title("Uji konvergensi skema Godunov (kondisi awal mulus)")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    path = os.path.join(OUT, "02_convergence.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[Konvergensi] Plot disimpan: {path}")
    return Ms, errors_l1, orders


# --------------------------------------------------------------------------
# Eksperimen 3: Kuantifikasi ketidakpastian Monte-Carlo
# --------------------------------------------------------------------------
def run_montecarlo(N: int = 400, M: int = 200, t_final: float = 1.0, seed: int = 42):
    """Bungkus solver dengan Monte-Carlo:
    - laju kedatangan (inflow density) acak
    - kapasitas mu acak (fluktuasi mesin)
    - amati throughput akhir (integral fluks di outflow) & lokasi kepadatan puncak.
    """
    rng = np.random.default_rng(seed)
    dx = 1.0 / M
    x = (np.arange(M) + 0.5) * dx

    throughput = np.empty(N)
    peak_x = np.empty(N)

    for k in range(N):
        # sampling stokastik
        lam = rng.poisson(lam=8) / 4.0            # laju kedatangan (0..~5)
        inflow_rho = 0.2 + 0.05 * lam             # konversi ke kerapatan hulu
        mu = float(rng.normal(loc=1.0, scale=0.15))
        mu = max(mu, 0.2)                          # jaga positif

        rho0 = np.full(M, 0.3)                    # WIP awal seragam ringan
        # jalankan solver deterministik
        rho = rho0.copy()
        t = 0.0
        flux_out_accum = 0.0
        while t < t_final:
            dt = cfl_dt(rho, dx, mu=mu, cfl=0.8)
            if t + dt > t_final:
                dt = t_final - t
            # fluks keluar di outflow = F(rho[-1])
            flux_out_accum += clearing_function(np.array([rho[-1]]), mu)[0] * dt
            rho = step_euler(rho, dx, dt, mu=mu, inflow=inflow_rho)
            t += dt

        throughput[k] = flux_out_accum
        peak_x[k] = x[int(np.argmax(rho))]

    # statistik
    mean_th, std_th = throughput.mean(), throughput.std(ddof=1)
    ci = 1.96 * std_th / np.sqrt(N)
    print(f"[Monte-Carlo] N={N}")
    print(f"[Monte-Carlo] Throughput: mean={mean_th:.4f}  std={std_th:.4f}  95% CI +-{ci:.4f}")
    print(f"[Monte-Carlo] Lokasi puncak WIP: mean={peak_x.mean():.3f}  std={peak_x.std(ddof=1):.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].hist(throughput, bins=30, color="steelblue", edgecolor="white")
    axes[0].axvline(mean_th, color="red", ls="--", label=f"mean = {mean_th:.3f}")
    axes[0].set_xlabel("throughput terintegrasi (unit)")
    axes[0].set_ylabel("frekuensi")
    axes[0].set_title(f"Distribusi throughput ({N} realisasi Monte-Carlo)")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].hist(peak_x, bins=30, color="darkorange", edgecolor="white")
    axes[1].set_xlabel("x lokasi WIP puncak (derajat penyelesaian)")
    axes[1].set_ylabel("frekuensi")
    axes[1].set_title("Distribusi lokasi puncak WIP")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    path = os.path.join(OUT, "03_montecarlo.png")
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"[Monte-Carlo] Plot disimpan: {path}")
    return throughput, peak_x


# --------------------------------------------------------------------------
# Eksperimen 4 (sanity): konservasi massa untuk sistem tertutup
# --------------------------------------------------------------------------
def run_mass_conservation(M: int = 200, t_final: float = 0.5, mu: float = 1.0):
    """Uji sederhana: tanpa inflow/outflow dinamis, cek massa total awal vs akhir.

    Karena batas hilir bersifat outflow bebas, kita justru mengharapkan massa
    berkurang seiring produk selesai. Test ini menampilkan neraca (massa awal,
    massa akhir, fluks keluar terintegrasi) -- keseimbangan harus terpenuhi:
        massa_akhir = massa_awal + inflow_terintegrasi - outflow_terintegrasi
    """
    dx = 1.0 / M
    x = (np.arange(M) + 0.5) * dx
    rho = np.exp(-((x - 0.5) ** 2) / 0.02)  # pulsa awal
    inflow_rho = 0.0  # tak ada kedatangan baru
    m0 = total_mass(rho, dx)

    t = 0.0
    outflow = 0.0
    while t < t_final:
        dt = cfl_dt(rho, dx, mu=mu, cfl=0.8)
        if t + dt > t_final:
            dt = t_final - t
        outflow += clearing_function(np.array([rho[-1]]), mu)[0] * dt
        rho = step_euler(rho, dx, dt, mu=mu, inflow=inflow_rho)
        t += dt
    m1 = total_mass(rho, dx)

    residu = m0 - m1 - outflow
    print(f"[Konservasi] massa awal   = {m0:.6f}")
    print(f"[Konservasi] massa akhir  = {m1:.6f}")
    print(f"[Konservasi] outflow int. = {outflow:.6f}")
    print(f"[Konservasi] residu (harus ~0) = {residu:+.2e}")
    return m0, m1, outflow, residu


if __name__ == "__main__":
    print("=" * 60)
    print("EKSPERIMEN 1: RIEMANN PROBLEM (SHOCK & RAREFACTION)")
    print("=" * 60)
    run_riemann()
    print()
    print("=" * 60)
    print("EKSPERIMEN 2: UJI KONVERGENSI")
    print("=" * 60)
    run_convergence()
    print()
    print("=" * 60)
    print("EKSPERIMEN 3: MONTE-CARLO (UQ)")
    print("=" * 60)
    run_montecarlo()
    print()
    print("=" * 60)
    print("EKSPERIMEN 4: NERACA MASSA")
    print("=" * 60)
    run_mass_conservation()
    print()
    print("Semua eksperimen selesai. Cek folder results/ untuk plot-nya.")
