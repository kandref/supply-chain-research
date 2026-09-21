"""
Deskripsi Lagrangian untuk model aliran produksi ADR.

Alih-alih menghitung rata-rata sel pada grid tetap (Eulerian), pendekatan
Lagrangian melacak sejumlah PAKET (parcel) yang masing-masing membawa massa
WIP tetap. Setiap paket k bergerak sepanjang sumbu penyelesaian menurut

    dX_k/dt = v(rho(X_k, t)),   v = F(rho)/rho = mu / (1 + rho)

Kerapatan makroskopik rho direkonstruksi dari sebaran posisi paket melalui
kernel density estimate (KDE) Gaussian, sehingga besaran Lagrangian dapat
dipetakan kembali ke grid Eulerian untuk perbandingan / coupling.
"""

from __future__ import annotations

import numpy as np

from adr_model import clearing_function


def init_parcels_from_density(rho0: np.ndarray, dx: float, parcels_per_unit_mass: int = 2000):
    """Inisialisasi posisi paket agar sebarannya mereproduksi rho0.

    Massa total = sum(rho0)*dx. Jumlah paket ~ parcels_per_unit_mass * massa.
    Tiap paket membawa massa sama (w). Posisi disebar proporsional terhadap
    rho0 memakai inversi CDF diskret.
    """
    M = rho0.size
    x_centers = (np.arange(M) + 0.5) * dx
    total_mass = float(np.sum(rho0) * dx)

    n_parcels = max(int(parcels_per_unit_mass * total_mass), 100)
    w = total_mass / n_parcels

    # CDF diskret dari rho0
    pmf = rho0 / np.sum(rho0)
    cdf = np.cumsum(pmf)

    # sampling posisi via inversi CDF + jitter dalam sel
    u = (np.arange(n_parcels) + 0.5) / n_parcels
    cell_idx = np.searchsorted(cdf, u)
    cell_idx = np.clip(cell_idx, 0, M - 1)
    rng = np.random.default_rng(0)
    jitter = (rng.random(n_parcels) - 0.5) * dx
    positions = x_centers[cell_idx] + jitter
    positions = np.clip(positions, 0.0, 1.0)

    return positions, w


def reconstruct_density(
    positions: np.ndarray,
    w: float,
    grid_x: np.ndarray,
    h: float,
    domain: tuple[float, float] = (0.0, 1.0),
):
    """Rekonstruksi rho pada grid dari sebaran paket via KDE Gaussian
    ter-renormalisasi di batas.

    rho(x_i) = [ sum_k w * K_h(x_i - x_k) ] / Z(x_i),
    K_h(r) = exp(-r^2 / (2 h^2)) / (sqrt(2 pi) h),
    Z(x_i) = integral_{a}^{b} K_h(x_i - x') dx'  (koreksi kebocoran di tepi).

    `domain=(a,b)` menyatakan batas fisik tempat massa dikekalkan. Faktor Z
    mengoreksi kernel yang "bocor" keluar [a,b] sehingga rho tidak terestimasi
    terlalu rendah di dekat tepi. Untuk coupling sub-window, domain diisi batas
    sub-window agar normalisasi benar (BUKAN selalu [0,1]).
    """
    if positions.size == 0:
        return np.zeros_like(grid_x)
    from math import erf

    a, b = domain
    diff = grid_x[:, None] - positions[None, :]
    kern = np.exp(-0.5 * (diff / h) ** 2) / (np.sqrt(2 * np.pi) * h)
    rho_raw = w * kern.sum(axis=1)

    s = np.sqrt(2.0) * h
    erf_vec = np.vectorize(erf)
    Z = 0.5 * (erf_vec((b - grid_x) / s) - erf_vec((a - grid_x) / s))
    Z = np.clip(Z, 1e-6, None)
    return rho_raw / Z


def interp_density_at(positions: np.ndarray, grid_x: np.ndarray, rho_grid: np.ndarray):
    """Interpolasi linear rho dari grid ke posisi paket."""
    return np.interp(positions, grid_x, rho_grid)


def solve_lagrangian(
    rho0: np.ndarray,
    t_final: float,
    dx: float,
    mu: float = 1.0,
    cfl: float = 0.4,
    parcels_per_unit_mass: int = 4000,
    bandwidth_factor: float = 2.0,
    inflow: float | None = None,
):
    """Integrasi model ADR secara Lagrangian.

    Mengembalikan (grid_x, rho_final) hasil rekonstruksi pada grid Eulerian
    yang sama, agar bisa dibandingkan langsung dengan solver Eulerian.
    """
    M = rho0.size
    grid_x = (np.arange(M) + 0.5) * dx
    h = bandwidth_factor * dx

    positions, w = init_parcels_from_density(rho0, dx, parcels_per_unit_mass)

    # kerapatan inflow di hulu (x=0): jika None, samakan dengan rho0[0]
    rho_inflow = float(rho0[0]) if inflow is None else float(inflow)

    t = 0.0
    inflow_residual = 0.0  # akumulator massa masuk pecahan (< w)
    while t < t_final:
        # rekonstruksi rho pada grid lalu interpolasi ke paket
        rho_grid = reconstruct_density(positions, w, grid_x, h)
        rho_at_parcels = interp_density_at(positions, grid_x, rho_grid)
        vel = mu / (1.0 + rho_at_parcels)

        # langkah waktu tunduk CFL berbasis kecepatan maksimum
        v_max = max(float(np.max(vel)), 1e-12) if vel.size else mu
        dt = cfl * dx / v_max
        if t + dt > t_final:
            dt = t_final - t

        positions = positions + vel * dt if positions.size else positions
        # outflow bebas: buang paket yang melewati x=1
        positions = positions[positions <= 1.0]

        # injeksi paket di hulu: massa masuk = F(rho_inflow) * dt
        inflow_residual += clearing_function(np.array([rho_inflow]), mu)[0] * dt
        n_new = int(inflow_residual / w)
        if n_new > 0:
            inflow_residual -= n_new * w
            # sebar paket baru di lapisan tipis dekat x=0 (lebar ~ h)
            new_pos = np.random.default_rng(int(t * 1e6) % (2**32)).random(n_new) * h
            positions = np.concatenate([positions, new_pos])

        t += dt

    rho_final = reconstruct_density(positions, w, grid_x, h)
    return grid_x, rho_final
