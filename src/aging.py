"""
Pelacakan umur (age) dan lead time per batch dalam model ADR.

Ini adalah peran BENAR untuk deskripsi Lagrangian dalam coupling multi-skala:
Eulerian menangani dinamika kerapatan makro (efisien & monoton di shock),
sementara Lagrangian berjalan paralel untuk melacak identitas dan riwayat
setiap batch produk -- besaran yang secara fundamental TIDAK DAPAT diperoleh
dari deskripsi Eulerian.

Yang dilacak per paket:
  - position (X_k)   : posisi penyelesaian saat ini
  - birth_time (t0_k): waktu paket masuk sistem (di x=0)

Dari data ini dapat dihitung:
  - age(x_i, t)      : distribusi umur WIP pada tahap x_i
  - lead_time (LT)   : t_keluar - t_masuk untuk paket yang mencapai x=1

Kecepatan paket ditentukan oleh kerapatan LOKAL (dari solver Eulerian yang
berjalan paralel) -> ini yang membuatnya "coupled".
"""

from __future__ import annotations

import numpy as np

from adr_model import (
    clearing_function,
    step_euler,
    cfl_dt,
)


def solve_aging_coupled(
    rho0: np.ndarray,
    t_final: float,
    dx: float,
    mu: float = 1.0,
    cfl: float = 0.5,
    inflow: float = 0.5,
    injection_rate: float | None = None,
    seed: int = 0,
):
    """Solver coupled Euler-Lagrange komplementer.

    - Eulerian menghitung rho(x, t) (kerapatan WIP)
    - Lagrangian: paket dengan (X, t_lahir) mengalir dengan kecepatan
      v = mu / (1 + rho_lokal), diinjeksi di hulu dengan laju yang sesuai
      dengan fluks masuk F(rho_inflow), dan dicabut saat X > 1.

    Return:
      grid_x, times, rho_hist, lead_times, ages_final, ages_at_all_x
    """
    rng = np.random.default_rng(seed)
    M = rho0.size
    grid_x = (np.arange(M) + 0.5) * dx

    rho = rho0.astype(float).copy()

    # inisialisasi paket agar sebaran awal ~ rho0 (representasi 1 paket = mass_per_parcel)
    mass_per_parcel = 0.01  # granularitas: 1 paket mewakili 0.01 unit massa
    m0 = float(np.sum(rho0) * dx)
    n0 = max(int(m0 / mass_per_parcel), 50)
    if np.sum(rho0) > 0:
        pmf = rho0 / np.sum(rho0)
        cdf = np.cumsum(pmf)
        u = (np.arange(n0) + 0.5) / n0
        idx = np.clip(np.searchsorted(cdf, u), 0, M - 1)
        positions = grid_x[idx] + (rng.random(n0) - 0.5) * dx
        positions = np.clip(positions, 0.0, 1.0)
    else:
        positions = np.array([])
    # Umur awal: paket dianggap uniformly-aged sesuai kerapatan setempat.
    # Untuk sederhana: t_lahir = 0 (semua paket awal "baru").
    birth_times = np.zeros_like(positions)
    # Tandai paket: True = diinjeksi dari hulu (x=0), False = seed awal.
    # Hanya paket ter-injeksi yang lead time-nya SAH (masuk penuh dari x=0);
    # paket seed mulai di tengah domain sehingga transit-nya bukan lead time.
    is_injected = np.zeros(positions.size, dtype=bool)

    lead_times = []  # akan diisi saat paket keluar di x >= 1

    rho_hist = [rho.copy()]
    times = [0.0]

    # laju injeksi paket di hulu = F(inflow) / mass_per_parcel [paket/waktu]
    if injection_rate is None:
        injection_rate = clearing_function(np.array([inflow]), mu)[0] / mass_per_parcel
    inject_residual = 0.0

    t = 0.0
    step = 0
    while t < t_final:
        dt = cfl_dt(rho, dx, mu=mu, cfl=cfl)
        if t + dt > t_final:
            dt = t_final - t

        # 1. update Eulerian (rho) satu langkah
        rho = step_euler(rho, dx, dt, mu=mu, inflow=inflow)

        # 2. update posisi paket: v = mu / (1 + rho_lokal)
        if positions.size > 0:
            rho_at_p = np.interp(positions, grid_x, rho)
            vel = mu / (1.0 + rho_at_p)
            positions = positions + vel * dt

            # paket yang keluar (X > 1) -> catat lead time HANYA jika ter-injeksi
            done_mask = positions > 1.0
            if np.any(done_mask):
                done_injected = done_mask & is_injected
                if np.any(done_injected):
                    lts = (t + dt) - birth_times[done_injected]
                    lead_times.extend(lts.tolist())
                positions = positions[~done_mask]
                birth_times = birth_times[~done_mask]
                is_injected = is_injected[~done_mask]

        # 3. injeksi paket baru di hulu
        inject_residual += injection_rate * dt
        n_new = int(inject_residual)
        if n_new > 0:
            inject_residual -= n_new
            # sebar tipis di lapisan dekat x=0 (setebal ~1 sel), waktu lahir = t
            new_pos = rng.random(n_new) * dx
            new_birth = np.full(n_new, t + 0.5 * dt)
            positions = np.concatenate([positions, new_pos])
            birth_times = np.concatenate([birth_times, new_birth])
            is_injected = np.concatenate([is_injected, np.ones(n_new, dtype=bool)])

        t += dt
        step += 1
        if step % 20 == 0:
            rho_hist.append(rho.copy())
            times.append(t)

    rho_hist.append(rho.copy())
    times.append(t)

    # snapshot umur akhir per paket
    ages_final = t - birth_times

    # bin umur ke grid untuk profil age(x) di akhir
    ages_by_x = np.full(M, np.nan)
    counts = np.zeros(M, dtype=int)
    sum_ages = np.zeros(M)
    if positions.size > 0:
        bin_idx = np.clip((positions / dx).astype(int), 0, M - 1)
        np.add.at(sum_ages, bin_idx, ages_final)
        np.add.at(counts, bin_idx, 1)
    nz = counts > 0
    ages_by_x[nz] = sum_ages[nz] / counts[nz]

    return (
        grid_x,
        np.array(times),
        np.array(rho_hist),
        np.array(lead_times),
        ages_final,
        ages_by_x,
        positions,
        birth_times,
    )


def littles_law_check(lead_times: np.ndarray, throughput_rate: float) -> tuple[float, float]:
    """Uji hukum Little: rata-rata WIP N = throughput * mean_lead_time.

    Return (N_prediksi, mean_lead_time). Nilai N_prediksi harus mendekati
    massa rata-rata di sistem pada steady state.
    """
    if lead_times.size == 0:
        return 0.0, 0.0
    L = float(np.mean(lead_times))
    N_pred = throughput_rate * L
    return N_pred, L
