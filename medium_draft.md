# Ketika Macet di Jalan Tol Mengajari Kita soal Pabrik: Memodelkan Bottleneck Rantai Pasok dengan Persamaan Gelombang

*Sebuah studi kasus numerik dengan Python — dari gelombang kejut hukum kekekalan sampai hukum Little.*

---

## Bermula dari sebuah kemacetan

Pernah memperhatikan bahwa kemacetan di jalan tol sering kali bergerak *mundur*? Mobil di depan sudah melaju, tapi gelombang rem-menyala terus merambat ke belakang, menabrak mobil-mobil yang baru datang. Anehnya, tidak ada kecelakaan, tidak ada penyebab yang jelas di titik macet itu sendiri. Kemacetan itu seolah hidup sendiri.

Fenomena ini punya nama: **gelombang kejut lalu lintas** (*traffic shock wave*), dan sudah lama dimodelkan dengan matematika yang sama yang dipakai untuk menggambarkan aliran fluida dan ledakan sonik — **hukum kekekalan hiperbolik** (*hyperbolic conservation laws*).

Yang menarik: struktur matematis yang sama ternyata muncul di tempat yang sama sekali berbeda — **lantai pabrik dan rantai pasok**. Ganti "mobil" dengan "produk setengah jadi", ganti "posisi di jalan" dengan "derajat penyelesaian produksi", dan tiba-tiba kemacetan tol berubah menjadi **bottleneck produksi**: tumpukan barang dalam proses (*work-in-process*, WIP) yang menggunung di depan sebuah stasiun kerja yang kapasitasnya terbatas.

Artikel ini menceritakan eksperimen numerik kecil yang saya bangun untuk mempelajari fenomena ini — lengkap dengan kode, verifikasi, dan beberapa temuan yang sempat membuat saya harus jujur pada diri sendiri.

---

## Model: kepadatan barang sebagai "fluida"

Alih-alih melacak tiap produk satu per satu (yang akan mahal secara komputasi untuk jutaan unit), kita perlakukan aliran produksi seperti fluida yang punya **kerapatan** ρ(x, t): berapa banyak WIP yang ada pada tahap penyelesaian x pada waktu t. Sumbu x di sini bukan posisi fisik, melainkan **derajat penyelesaian** dari 0 (bahan mentah) sampai 1 (produk jadi).

Dinamikanya diatur persamaan kekekalan:

> ∂ρ/∂t + ∂F(ρ)/∂x = 0

dengan F(ρ) adalah **clearing function** — laju keluaran produksi sebagai fungsi seberapa padat WIP. Saya memakai bentuk yang populer di literatur rantai pasok (model Armbruster–Degond–Ringhofer):

> F(ρ) = μ·ρ / (1 + ρ)

Interpretasinya intuitif: saat WIP sedikit, output naik sebanding beban; saat WIP menumpuk, output jenuh ke kapasitas maksimum μ. Ini persis perilaku mesin nyata.

---

## Langkah 1: Membuktikan solvernya tidak bohong

Sebelum percaya pada hasil apa pun, solver numerik harus diverifikasi terhadap kasus yang jawabannya sudah kita tahu. Ini bagian yang sering dilewati orang, padahal paling penting.

**Riemann problem** — dua kepadatan konstan bertemu di satu titik — adalah batu uji klasik. Teori (kondisi Rankine–Hugoniot) memberi tahu kita persis seberapa cepat gelombang kejut harus merambat. Solver saya (skema *finite volume* Godunov) menghasilkan posisi shock yang meleset hanya **6×10⁻⁴** dari nilai analitis.

![Riemann problem: pembentukan gelombang kejut dan penjarangan](IMG_01_riemann)

Diagram ruang-waktu memperlihatkan "kemacetan" itu merambat sebagai garis diagonal yang tajam — persis analog dengan gelombang rem di tol tadi.

![Diagram ruang-waktu: bottleneck merambat dengan kecepatan Rankine–Hugoniot](IMG_01_spacetime)

Uji kedua: **konvergensi**. Kalau grid diperhalus, galat harus mengecil dengan laju yang bisa diprediksi. Hasilnya orde ≈ 1.0, persis sesuai teori skema Godunov orde satu.

![Uji konvergensi: galat L1 mengecil sesuai orde teoretis](IMG_02_convergence)

Dan uji **konservasi massa**: barang tidak boleh hilang atau muncul dari ketiadaan. Residu neraca massa saya **~7×10⁻¹⁸** — setara batas presisi komputer. Solver lolos.

---

## Langkah 2: Dunia itu tidak pasti — Monte-Carlo

Pabrik nyata tidak deterministik. Kedatangan order berfluktuasi (distribusi Poisson), kapasitas mesin bergoyang. Dengan menjalankan **400 realisasi** dengan parameter acak, kita dapat bukan satu angka, tapi seluruh *distribusi* kemungkinan throughput — informasi yang jauh lebih berguna untuk mengambil keputusan.

![Distribusi throughput dari 400 simulasi Monte-Carlo](IMG_03_montecarlo)

---

## Langkah 3: Dua sudut pandang, dan sebuah pelajaran pahit

Ada dua cara memandang aliran: **Eulerian** (kita duduk diam, mengamati kepadatan lewat) dan **Lagrangian** (kita ikut menaiki tiap paket produk). Saya membangun keduanya dan mencocokkannya. Pada kondisi mulus, keduanya sepakat sampai galat 2×10⁻².

![Perbandingan silang Eulerian vs Lagrangian pada kondisi mulus](IMG_04_xcheck)

Lalu saya berasumsi — dengan percaya diri — bahwa deskripsi Lagrangian akan **lebih tajam** di sekitar bottleneck. Ternyata **saya salah**.

![Coupling di shock: Godunov justru lebih tajam, Lagrangian overshoot](IMG_05_coupling)

Skema Godunov menghasilkan shock dengan lebar transisi 0.005, sementara rekonstruksi Lagrangian malah *overshoot* sampai 3.8 dengan osilasi. Ini bukan bug yang bisa ditambal — ini sifat fundamental: skema Godunov sudah monoton (TVD) dan memang unggul di diskontinuitas.

Alih-alih memaksakan asumsi awal, saya merumuskan ulang perannya. **Lagrangian tidak untuk menajamkan shock. Lagrangian untuk melacak sesuatu yang Eulerian secara struktural tidak bisa: identitas dan riwayat tiap batch.**

---

## Langkah 4: Umur dan lead time — di sinilah Lagrangian bersinar

Deskripsi Eulerian tahu *berapa banyak* WIP di tiap tahap, tapi tidak tahu **sudah berapa lama** sebuah batch tertentu berada di sistem. Padahal itulah metrik yang dipedulikan manajer pabrik: **lead time**.

Dengan menandai tiap paket dengan waktu lahirnya, saya bisa merekam distribusi lead time produk yang selesai. Hasilnya: rata-rata lead time **1.4914** vs prediksi analitis **1.5000** — meleset 0.6%.

Uji paling meyakinkan adalah **hukum Little** (N = λ · LT), hukum universal antrian. Rasio prediksi terhadap WIP aktual: **0.994**. Karena hukum Little berlaku untuk sistem tunak apa pun, kesesuaian ini adalah bukti kuat bahwa pelacakannya benar secara fisis.

![Pelacakan umur & lead time: profil WIP, distribusi lead time, profil umur](IMG_06_aging)

---

## Langkah 5: Menciptakan bottleneck, dan menyaksikan kekacauan

Sekarang bagian yang paling seru. Saya turunkan kapasitas satu tahap ke 35% pada jendela x ∈ [0.55, 0.65], meniru sebuah stasiun kerja yang kewalahan. Untuk menghasilkan *back-pressure* (antrian merambat mundur, persis gelombang tol), saya perluas flux ke bentuk gaya **Cell Transmission Model**.

Hasilnya dramatis:

- **Lead time rata-rata**: naik **4×** (1.98 → 7.91)
- **Lead time persentil-95**: naik **6.6×** (2.0 → 13.13)
- **Variabilitas lead time**: naik **46×**
- **WIP total**: naik **3.7×**
- **Throughput**: turun ~51%

![Dampak bottleneck: WIP menumpuk, lead time berekor panjang, umur melonjak](IMG_07_bottleneck)

Perhatikan tiga hal di gambar: (1) WIP menumpuk jenuh di *hulu* bottleneck — buffer penuh karena barang tak bisa lewat; (2) distribusi lead time berubah dari paku tajam menjadi ekor panjang penuh ketidakpastian; (3) umur WIP melonjak tepat di depan bottleneck. Ini adalah cerita kuantitatif tentang mengapa satu stasiun yang lambat bisa melumpuhkan seluruh pabrik.

---

## Penutup: matematika yang sama, dunia yang berbeda

Yang saya sukai dari proyek ini adalah betapa satu kerangka matematika — hukum kekekalan hiperbolik yang lahir dari dinamika gas dan lalu lintas — bisa berbicara dengan fasih tentang lantai pabrik. Dan bagaimana dua sudut pandang (Eulerian untuk kepadatan makro, Lagrangian untuk riwayat individual) ternyata **saling melengkapi**, bukan bersaing.

Pelajaran terbesar buat saya justru bukan soal matematikanya, tapi soal kejujuran: asumsi awal saya soal Lagrangian ternyata keliru, dan mengakuinya malah membuka arah kontribusi yang jauh lebih menarik.

Kode lengkap, semua eksperimen, dan proposal risetnya terbuka di:
**https://github.com/kandref/supply-chain-research**

Kalau kamu tertarik ngobrol soal pemodelan numerik, rantai pasok, atau kenapa kemacetan tol itu indah secara matematis — sila sapa. 👋

---

*Ditulis oleh Kurnia Andre Febrian. Dibangun dengan Python, NumPy, dan Matplotlib.*
