"""
Model aliran produksi Armbruster-Degond-Ringhofer (ADR).

Menyelesaikan hukum kekekalan hiperbolik 1D:

    d(rho)/dt + d(F(rho))/dx = 0,   x in [0, 1]

dengan x = derajat penyelesaian produk (0 = bahan mentah, 1 = barang jadi)
dan rho = kerapatan barang dalam proses (WIP).

Clearing function (analog fundamental diagram lalu lintas):

    F(rho) = mu * rho / (1 + rho)

Skema: finite volume + fluks Godunov, langkah waktu adaptif tunduk CFL.
"""

from __future__ import annotations

import numpy as np


def clearing_function(rho: np.ndarray, mu: float = 1.0) -> np.ndarray:
    """Clearing function F(rho) = mu * rho / (1 + rho).

    Cekung (concave), monoton naik, jenuh menuju mu saat rho -> inf.
    """
    return mu * rho / (1.0 + rho)


def flux_derivative(rho: np.ndarray, mu: float = 1.0) -> np.ndarray:
    """Turunan F'(rho) = mu / (1 + rho)^2 -- kecepatan karakteristik."""
    return mu / (1.0 + rho) ** 2


def godunov_flux(rho_l: np.ndarray, rho_r: np.ndarray, mu: float = 1.0) -> np.ndarray:
    """Fluks numerik Godunov pada batas antar-sel.

    Untuk clearing function yang CEKUNG dan MONOTON NAIK, F'(rho) > 0 untuk
    semua rho >= 0, sehingga tidak ada titik kritis interior: kecepatan
    karakteristik selalu positif (arah aliran hulu -> hilir). Solusi Riemann
    Godunov untuk flux cekung:

        rho_l <= rho_r : F_hat = min(F(rho_l), F(rho_r))
        rho_l >  rho_r : F_hat = max_{rho in [rho_r, rho_l]} F(rho)

    Karena F monoton naik, min = F(rho_l) dan max = F(rho_l); jadi untuk kasus
    ini fluks Godunov tereduksi menjadi upwind F(rho_l). Kita tetap tulis bentuk
    umumnya agar kokoh bila clearing function diganti yang punya titik kritis.
    """
    fl = clearing_function(rho_l, mu)
    fr = clearing_function(rho_r, mu)

    # F monoton naik pada [0, inf) => tidak ada rho kritis interior.
    # min pada [rho_l, rho_r] ada di ujung kiri; max pada [rho_r, rho_l] di ujung kiri.
    f_upwind = np.where(rho_l <= rho_r, np.minimum(fl, fr), np.maximum(fl, fr))
    return f_upwind


def step_euler(
    rho: np.ndarray,
    dx: float,
    dt: float,
    mu: float = 1.0,
    inflow: float | None = None,
) -> np.ndarray:
    """Satu langkah waktu skema finite volume (deskripsi Eulerian).

    rho     : rata-rata sel saat ini (shape [M])
    inflow  : kerapatan yang masuk di batas hulu (x=0). Jika None, pakai rho[0].
    """
    m = rho.size
    # bentuk state dengan ghost cell di kedua ujung
    rho_l = np.empty(m + 1)  # nilai kiri tiap batas (M+1 batas)
    rho_r = np.empty(m + 1)  # nilai kanan tiap batas

    # batas interior i+1/2: kiri = sel i, kanan = sel i+1
    rho_l[1:] = rho          # batas 1..M kiri = sel 0..M-1
    rho_r[:-1] = rho         # batas 0..M-1 kanan = sel 0..M-1

    # batas hulu (index 0): kiri = inflow (ghost), kanan = sel 0
    rho_l[0] = inflow if inflow is not None else rho[0]
    # batas hilir (index M): kanan = sel M-1 (outflow bebas, ghost = sel terakhir)
    rho_r[-1] = rho[-1]

    fluxes = godunov_flux(rho_l, rho_r, mu)  # shape [M+1]

    rho_new = rho - (dt / dx) * (fluxes[1:] - fluxes[:-1])
    # jaga non-negatif secara numerik
    return np.maximum(rho_new, 0.0)


def cfl_dt(rho: np.ndarray, dx: float, mu: float = 1.0, cfl: float = 0.8) -> float:
    """Langkah waktu adaptif dari syarat CFL: dt = cfl * dx / max|F'(rho)|."""
    a_max = np.max(np.abs(flux_derivative(rho, mu)))
    a_max = max(a_max, 1e-12)
    return cfl * dx / a_max


def solve_adr(
    rho0: np.ndarray,
    t_final: float,
    dx: float,
    mu: float = 1.0,
    cfl: float = 0.8,
    inflow: float | None = None,
    save_every: int = 1,
):
    """Integrasi model ADR dari t=0 sampai t_final.

    Mengembalikan (times, history) dengan history[k] = profil rho pada times[k].
    """
    rho = rho0.astype(float).copy()
    t = 0.0
    times = [0.0]
    history = [rho.copy()]
    step = 0

    while t < t_final:
        dt = cfl_dt(rho, dx, mu, cfl)
        if t + dt > t_final:
            dt = t_final - t
        rho = step_euler(rho, dx, dt, mu, inflow)
        t += dt
        step += 1
        if step % save_every == 0:
            times.append(t)
            history.append(rho.copy())

    if times[-1] < t_final - 1e-12:
        times.append(t)
        history.append(rho.copy())

    return np.array(times), np.array(history)


def total_mass(rho: np.ndarray, dx: float) -> float:
    """Massa total WIP (integral rho dx) -- harus lestari tanpa inflow/outflow."""
    return float(np.sum(rho) * dx)
