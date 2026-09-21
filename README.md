# supply-chain-research

Riset numerik untuk fenomena **bottleneck** pada aliran produksi rantai pasok
menggunakan model kekekalan hiperbolik (Armbruster–Degond–Ringhofer) dan
kuantifikasi ketidakpastian berbasis Monte-Carlo.

Repositori ini berisi:

- Proposal S2 dalam LaTeX (folder `proposal/`, cover → Bab 3).
- Solver numerik model ADR dalam Python (folder `src/`).
- Eksperimen verifikasi + plot hasil (folder `results/`).

## Model

Model ADR menyatakan aliran produksi sebagai hukum kekekalan hiperbolik satu
dimensi

```
d(rho)/dt + d(F(rho))/dx = 0,   x in [0, 1]
```

dengan `x` = derajat penyelesaian produk (0 = bahan mentah, 1 = barang jadi),
`rho` = kerapatan WIP (*work in progress*), dan *clearing function*

```
F(rho) = mu * rho / (1 + rho)
```

Skema numerik: **finite volume** + fluks **Godunov**, langkah waktu adaptif
CFL.

## Struktur

```
supply-chain-research/
├── src/
│   ├── adr_model.py       # solver deterministik (fluks Godunov, CFL)
│   └── experiments.py     # 4 eksperimen verifikasi
├── results/               # plot PNG hasil eksperimen
├── proposal/              # dokumen LaTeX (opsional, mengikuti PDF proposal)
├── requirements.txt
└── README.md
```

## Menjalankan

```bash
# 1. buat venv & install
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 2. jalankan semua eksperimen (menghasilkan plot di results/)
cd src
python experiments.py
```

Atau pakai `uv` yang lebih ringkas:

```bash
uv venv --python 3.11 .venv
source .venv/bin/activate
uv pip install -r requirements.txt
python src/experiments.py
```

## Eksperimen yang dijalankan

### 1. Riemann problem — pembentukan & pelepasan bottleneck

Kondisi awal berupa dua keadaan konstan yang dipisahkan diskontinuitas.
Skema harus mereproduksi:

- **Shock**: penumpukan WIP merambat mundur (bottleneck menjalar).
- **Rarefaction**: penurunan WIP menyebar mulus (bottleneck teratasi).

Kecepatan rambat shock dibandingkan dengan nilai analitis Rankine–Hugoniot

```
s = [F(rho_R) - F(rho_L)] / [rho_R - rho_L]
```

**Hasil**: selisih posisi shock numerik vs analitis pada t = 0.6 sebesar
**~6×10⁻⁴** (M = 400 sel).

### 2. Uji konvergensi

Kondisi awal mulus `rho(x,0) = 1 + 0.5 sin(2πx)`, dijalankan pada
M ∈ {50, 100, 200, 400}, dibandingkan dengan solusi referensi (M = 1600).

**Hasil**: orde konvergensi empiris `p ≈ 1.0`, sesuai orde teoretis skema
Godunov orde-1.

### 3. Monte-Carlo — kuantifikasi ketidakpastian

Membungkus solver dengan N realisasi acak:

- Laju kedatangan `λ ~ Poisson`
- Kapasitas `μ ~ Normal(1.0, 0.15)`

Keluaran: distribusi *throughput* terintegrasi dan lokasi puncak WIP.

**Hasil** (N = 400): throughput mean = 0.230, std = 0.037, 95% CI ±0.004.

### 4. Neraca massa

Uji sanity: massa awal − massa akhir − outflow terintegrasi ≈ 0 (sampai
presisi mesin). Verifikasi bahwa skema *finite volume* memang konservatif
secara diskret.

**Hasil**: residu **~7×10⁻¹⁸** (presisi mesin).

## Hasil visual

Semua plot ada di `results/`:

- `01_riemann.png` — snapshot profil WIP (shock & rarefaction)
- `01_riemann_spacetime.png` — diagram ruang-waktu shock
- `02_convergence.png` — grafik loglog galat vs Δx
- `03_montecarlo.png` — histogram throughput & lokasi puncak

## Roadmap

- [ ] Deskripsi Lagrangian (particle-based) + coupling multi-skala
- [ ] Multi-level Monte-Carlo untuk efisiensi UQ
- [ ] Validasi data produksi (opsional)
- [ ] Paralelisasi realisasi Monte-Carlo

## Lisensi

MIT.
