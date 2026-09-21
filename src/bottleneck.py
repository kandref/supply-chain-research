"""
Skenario bottleneck: kapasitas turun di satu tahap produksi.

CATATAN MODEL (penting & jujur):
Clearing function ADR F(rho) = mu*rho/(1+rho) bersifat MONOTON NAIK, sehingga
kecepatan karakteristik F'(rho) > 0 selalu -- informasi hanya merambat hilir.
Konsekuensinya, model ADR murni TIDAK dapat menghasilkan antrian (queue) yang
merambat MUNDUR ke hulu akibat bottleneck, dan penumpukan pada satu sel dapat
tumbuh tak terbatas tanpa "tekanan balik" (back-pressure).

Untuk menangkap back-pressure yang realistis dan menjaga solusi terbatas, di
sini flux diperluas ke bentuk gaya Cell Transmission Model (CTM):

    Demand  (sending)  D_i = mu_i * rho_i / (1 + rho_i)        (<= mu_i)
    Supply  (receiving) S_i = mu_i * clip(1 - rho_i/rho_max, 0, 1)
    Flux antarmuka   f_{i+1/2} = min( D_i, S_{i+1} )

Buffer kapasitas rho_max membatasi rho di tiap tahap; ketika sel penuh, supply
turun -> flux masuk berkurang -> antrian merambat mundur. Ini ekstensi baku
(Daganzo 1994) yang cocok untuk aliran produksi dengan buffer terbatas.

mu_i (kapasitas per tahap) dibuat turun di jendela bottleneck.
"""

from __future__ import annotations

import numpy as np


def capacity_field(M: int, dx: float, mu_normal: float, mu_bottleneck: float,
                   x_lo: float, x_hi: float) -> np.ndarray:
    """Profil kapasitas mu(x): mu_normal, turun ke mu_bottleneck pada [x_lo, x_hi]."""
    x = (np.arange(M) + 0.5) * dx
    mu = np.full(M, mu_normal)
    mask = (x >= x_lo) & (x <= x_hi)
    mu[mask] = mu_bottleneck
    return mu


def ctm_fluxes(rho, mu, rho_max, q_in):
    """Hitung flux antarmuka gaya CTM. Return array shape [M+1]."""
    D = mu * rho / (1.0 + rho)                      # demand / sending
    S = mu * np.clip(1.0 - rho / rho_max, 0.0, 1.0)  # supply / receiving
    M = rho.size
    f = np.empty(M + 1)
    f[0] = min(q_in, S[0])            # inflow boundary (dibatasi supply sel 0)
    f[1:M] = np.minimum(D[:-1], S[1:])  # antarmuka interior
    f[M] = D[-1]                       # outflow bebas (hilir kapasitas tak hingga)
    return f, D, S


def solve_bottleneck_aging(
    M: int = 200,
    t_final: float = 20.0,
    mu_normal: float = 1.0,
    mu_bottleneck: float = 0.35,
    x_lo: float = 0.55,
    x_hi: float = 0.65,
    q_in: float = 0.5,
    rho_max: float = 8.0,
    cfl: float = 0.4,
    mass_per_parcel: float = 0.01,
    seed: int = 0,
):
    """Solver Eulerian CTM + pelacakan Lagrangian umur/lead-time, dengan bottleneck.

    Return dict berisi grid, profil akhir, lead_times, ages_by_x, mu_field, dll.
    """
    rng = np.random.default_rng(seed)
    dx = 1.0 / M
    grid_x = (np.arange(M) + 0.5) * dx
    mu = capacity_field(M, dx, mu_normal, mu_bottleneck, x_lo, x_hi)

    rho = np.full(M, 0.02)  # mulai hampir kosong

    # partikel Lagrangian
    positions = np.array([])
    birth_times = np.array([])
    is_injected = np.array([], dtype=bool)
    lead_times = []

    inject_residual = 0.0
    t = 0.0
    step = 0
    rho_snaps = [rho.copy()]
    snap_times = [0.0]

    # kecepatan karakteristik maks (demand & supply branch)
    a_max = max(mu_normal, mu_normal / rho_max, mu_normal)  # ~ mu_normal
    dt = cfl * dx / a_max

    while t < t_final:
        if t + dt > t_final:
            dt = t_final - t

        f, D, S = ctm_fluxes(rho, mu, rho_max, q_in)
        rho = rho + (dt / dx) * (f[:-1] - f[1:])
        rho = np.clip(rho, 0.0, rho_max)

        # kecepatan transport per sel = flux keluar / rho (v = f_out / rho)
        f_out = f[1:]
        with np.errstate(divide="ignore", invalid="ignore"):
            v_cell = np.where(rho > 1e-9, f_out / rho, mu / (1.0 + rho))

        # update partikel
        if positions.size > 0:
            v_at_p = np.interp(positions, grid_x, v_cell)
            positions = positions + v_at_p * dt
            done = positions > 1.0
            if np.any(done):
                di = done & is_injected
                if np.any(di):
                    lead_times.extend(((t + dt) - birth_times[di]).tolist())
                positions = positions[~done]
                birth_times = birth_times[~done]
                is_injected = is_injected[~done]

        # injeksi paket di hulu sesuai flux masuk aktual f[0]
        inject_residual += f[0] * dt / mass_per_parcel
        n_new = int(inject_residual)
        if n_new > 0:
            inject_residual -= n_new
            new_pos = rng.random(n_new) * dx
            positions = np.concatenate([positions, new_pos])
            birth_times = np.concatenate([birth_times, np.full(n_new, t + 0.5 * dt)])
            is_injected = np.concatenate([is_injected, np.ones(n_new, dtype=bool)])

        t += dt
        step += 1
        if step % 40 == 0:
            rho_snaps.append(rho.copy())
            snap_times.append(t)

        # dt adaptif ringan: kecepatan bisa turun di queue, tapi a_max tetap aman
        dt = cfl * dx / a_max

    rho_snaps.append(rho.copy())
    snap_times.append(t)

    # profil umur akhir per x
    ages_final = t - birth_times
    ages_by_x = np.full(M, np.nan)
    if positions.size > 0:
        counts = np.zeros(M)
        sums = np.zeros(M)
        idx = np.clip((positions / dx).astype(int), 0, M - 1)
        np.add.at(sums, idx, ages_final)
        np.add.at(counts, idx, 1)
        nz = counts > 0
        ages_by_x[nz] = sums[nz] / counts[nz]

    return {
        "grid_x": grid_x,
        "mu": mu,
        "rho_final": rho,
        "rho_snaps": np.array(rho_snaps),
        "snap_times": np.array(snap_times),
        "lead_times": np.array(lead_times),
        "ages_by_x": ages_by_x,
        "dx": dx,
        "x_lo": x_lo,
        "x_hi": x_hi,
    }
