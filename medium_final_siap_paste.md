# Ketika Macet di Jalan Tol Mengajari Kita Soal Pabrik

_Cerita tentang memodelkan bottleneck rantai pasok dengan persamaan gelombang, dan tentang menghidupkan kembali sebuah proposal yang lima tahun terkubur._

## Sedikit cerita, sebelum masuk ke matematika

Tahun 2021 aku sempat hampir mendaftar S2 Fisika Teoretis di salah satu kampus, entah kampus S1-ku dulu atau kampus yang berkali-kali aku incar saat S1 namun urung didapatkan. Judul proposalnya panjang dan agak sok penting, kira-kira "Analisis Penggunaan CFD untuk Pemetaan dan Rekayasa Lalu Lintas." Aku bikin sampai cover, lalu keburu dapat kerja dan proposalnya masuk laci. Selesai.

Lima tahun kemudian, di tengah kerjaan sehari-hari sebagai analis BI, aku iseng buka file lama itu. Tiba-tiba pengin coba lagi, tapi dengan sudut pandang yang lebih dekat ke kerjaan sekarang: rantai pasok, bukan lalu lintas. Bottleneck produksi, bukan macet tol. Meskipun secara matematika, jujurnya, dua hal itu adalah adik-kakak.

Yang menarik, riset ulang ini aku kerjakan bareng sebuah agen AI bernama Hermes Agent, agent framework buatan Nous Research yang bisa nulis kode, jalanin eksperimen numerik, compile LaTeX, dan pada saat yang sama juga jujur waktu asumsiku ternyata keliru (dan aku akan cerita soal itu di bawah). Bukan sekadar auto-complete. Lebih mirip teman ngoprek yang tahan begadang: aku dikte arah dan justifikasi ilmiahnya, dia bantu turunkan jadi kode, plot, dan tulisan. Iterasi jadi cepat sekali, yang biasanya butuh berminggu-minggu buat setup solver dan verifikasi, selesai dalam beberapa sesi obrolan.

Aku tulis ini sebagian sebagai catatan teknis, sebagian sebagai catatan pribadi. Kalau kamu tertarik ke sisi teknisnya saja, lompati bagian ini dan langsung ke bawah.

## Kenapa pabrik itu mirip jalan tol

Pernah lihat macet di tol yang bergerak mundur? Mobil di depan sudah jalan, tapi gelombang rem terus merambat ke belakang, menabrak mobil yang baru datang. Fenomena ini punya nama, gelombang kejut lalu lintas (traffic shockwave), dan sudah lama dimodelkan dengan matematika yang sama yang dipakai buat aliran fluida.

Kejutannya, struktur matematis yang sama muncul juga di lantai pabrik. Ganti "mobil" dengan "produk setengah jadi," ganti "posisi di jalan" dengan "tahap penyelesaian produk," dan tiba-tiba kemacetan tol berubah jadi tumpukan barang di depan sebuah mesin yang kapasitasnya terbatas.

Alih-alih melacak jutaan unit produk satu per satu, kita perlakukan aliran produksi itu seperti fluida yang punya kerapatan ρ(x, t), berapa banyak barang dalam proses (WIP) yang ada pada tahap x pada waktu t. Sumbu x di sini bukan posisi fisik, melainkan "derajat penyelesaian," dari 0 (bahan mentah) sampai 1 (produk jadi).

Dinamikanya diatur satu persamaan sederhana, hukum kekekalan hiperbolik:

> ∂ρ/∂t + ∂F(ρ)/∂x = 0

F(ρ) adalah clearing function, laju keluaran produksi sebagai fungsi seberapa padat WIP. Aku pakai bentuk yang populer di literatur rantai pasok, model Armbruster, Degond, Ringhofer:

> F(ρ) = μ · ρ / (1 + ρ)

Kalau WIP sedikit, output naik sebanding beban. Kalau WIP menumpuk, output jenuh ke kapasitas maksimum μ. Ini persis perilaku mesin nyata yang kadang bisa cepat, kadang macet.

## Membuktikan solvernya tidak bohong

Sebelum percaya hasil apa pun, solver numerik harus dites. Ini bagian yang sering dilewati orang, padahal paling penting, aku sendiri nyaris melewatinya sebelum diingatkan.

Tiga uji yang aku jalankan:

- Riemann problem, batu uji klasik: dua kepadatan konstan bertemu di satu titik, dan teori (kondisi Rankine-Hugoniot) sudah kasih tahu persis seberapa cepat gelombang kejut harus merambat. Solverku (skema finite volume Godunov) menghasilkan posisi shock yang meleset 0.0006 dari nilai analitis. Kecil sekali.
- Uji konvergensi, kalau grid diperhalus, galat harus mengecil dengan laju yang bisa diprediksi. Hasilnya orde 1.0, persis sesuai teori skema Godunov orde satu.
- Konservasi massa, barang tidak boleh hilang atau muncul dari ketiadaan. Residu neraca massaku 7 × 10⁻¹⁸, setara batas presisi mesin.

Diagram ruang-waktu bikin fenomenanya lebih jelas. Gelombang kejut merambat sebagai garis diagonal yang tajam, persis analog dengan gelombang rem di tol tadi.


╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/01_riemann_spacetime.png        ║
║  Caption: Diagram ruang-waktu bottleneck yang merambat;    ║
║  garis diskontinuitas sesuai kecepatan Rankine-Hugoniot.   ║
╚══════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/01_riemann.png   (opsional)     ║
║  Caption: Riemann problem: pembentukan gelombang kejut     ║
║  (kiri) dan gelombang penjarangan (kanan).                 ║
╚══════════════════════════════════════════════════════════╝

╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/02_convergence.png  (opsional)  ║
║  Caption: Uji konvergensi: galat L1 mengecil terhadap Δx   ║
║  (skala log-log), orde empiris mendekati 1.0.              ║
╚══════════════════════════════════════════════════════════╝

Solver lolos ketiga ujian. Boleh lega, dan boleh lanjut.

## Dunia itu tidak pernah pasti

Pabrik nyata tidak deterministik. Order berfluktuasi, kapasitas mesin bergoyang, operator sakit, bahan baku telat. Kalau aku cuma jalankan simulasi sekali, hasilnya kelihatan meyakinkan padahal cuma satu skenario dari ribuan kemungkinan.

Jadi aku jalankan 400 kali dengan parameter acak, kedatangan order pakai distribusi Poisson, kapasitas pakai distribusi normal. Hasilnya bukan satu angka throughput, tapi seluruh distribusi. Ini yang bikin model jadi berguna buat pengambil keputusan, karena mereka bisa lihat bukan cuma "rata-rata," tapi juga "seberapa buruk kemungkinan terburuknya."


╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/03_montecarlo.png              ║
║  Caption: Distribusi throughput dari 400 realisasi        ║
║  Monte-Carlo, bukan satu angka tapi seluruh sebaran.      ║
╚══════════════════════════════════════════════════════════╝

## Dua sudut pandang, dan pengakuan bahwa aku salah

Ada dua cara memandang aliran. Eulerian, di mana kita duduk diam di satu titik dan mengamati kepadatan yang lewat. Lagrangian, di mana kita naik ke atas tiap paket produk dan ikut jalan bersamanya. Aku bangun keduanya dan cocokkan.

Pada kondisi mulus, keduanya sepakat sampai galat 2 persen. Bagus.


╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/04_xcheck_eul_lag.png          ║
║  Caption: Perbandingan silang Eulerian vs Lagrangian      ║
║  pada kondisi awal mulus; keduanya sepakat.               ║
╚══════════════════════════════════════════════════════════╝

Lalu aku berasumsi, dengan cukup percaya diri, bahwa deskripsi Lagrangian akan lebih tajam di sekitar bottleneck. Alasannya intuitif: kalau kita ikut naik ke tiap paket, kita tahu persis di mana tiap paket berada, jadi transisi tajam mestinya lebih mudah tertangkap.

Ternyata aku salah.


╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/05_coupling.png                ║
║  Caption: Di sekitar shock, Godunov justru lebih tajam;   ║
║  Lagrangian overshoot sampai 3.8 dengan osilasi.          ║
╚══════════════════════════════════════════════════════════╝

Skema Godunov menghasilkan shock dengan lebar transisi 0.005, sementara rekonstruksi Lagrangian malah overshoot sampai kerapatan 3.8, dengan osilasi. Ini bukan bug yang bisa ditambal, sifat teoretisnya memang begitu. Godunov punya properti monotonisitas yang bikin dia unggul di diskontinuitas, dan kernel rekonstruksi Lagrangian punya kelemahan alami di titik transisi tajam.

Alih-alih memaksakan asumsi lama, aku rumuskan ulang perannya. Lagrangian bukan untuk menajamkan shock. Lagrangian untuk melacak sesuatu yang Eulerian secara struktural tidak bisa: identitas dan riwayat tiap batch.

> Momen ini yang paling berkesan buat aku selama ngerjain proyek ini. Waktu Hermes bantu diagnosa hasil yang aneh itu, dia tidak berusaha buat asumsi awalku terlihat benar. Dia bilang, dengan bahasa yang lebih halus, bahwa premis awalku keliru, dan justru mengarahkan aku ke reframing yang jauh lebih menarik. Kejujuran seperti ini yang aku butuhkan.

## Umur dan lead time: di sinilah Lagrangian menemukan tempatnya

Deskripsi Eulerian tahu berapa banyak WIP di tiap tahap, tapi tidak tahu sudah berapa lama sebuah batch tertentu berada di sistem. Padahal itulah metrik yang manajer pabrik pedulikan setiap hari: lead time, umur produk, umur bahan baku sebelum terpakai.

Dengan menandai tiap paket dengan waktu lahirnya, aku bisa merekam distribusi lead time produk yang selesai. Rata-ratanya keluar di 1.4914, sedangkan prediksi analitisnya 1.5000, meleset 0.6 persen.

Tapi validasi yang paling meyakinkan datang dari hukum Little: jumlah WIP di sistem sama dengan laju throughput dikali rata-rata lead time. Hukum ini berlaku universal untuk sistem antrian tunak apa pun, dari kasir minimarket sampai kilang minyak. Rasio prediksi hukum Little terhadap WIP aktual di simulasiku, 0.994, meleset kurang dari 1 persen. Karena hukum ini begitu umum, kesesuaian ini adalah bukti kuat bahwa pelacakannya benar secara fisis.


╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/06_aging_leadtime.png          ║
║  Caption: Pelacakan umur dan lead time: profil WIP tunak  ║
║  (kiri), distribusi lead time (tengah), profil umur WIP   ║
║  sepanjang tahap penyelesaian (kanan).                    ║
╚══════════════════════════════════════════════════════════╝

## Menciptakan bottleneck dan menyaksikan kekacauan

Sekarang bagian paling seru. Aku turunkan kapasitas satu tahap ke 35 persen pada jendela sempit di x antara 0.55 sampai 0.65, simulasi dari satu stasiun kerja yang kewalahan. Buat menghasilkan antrian yang merambat mundur (persis fenomena rem tol tadi), aku perluas flux ke bentuk Cell Transmission Model, ekstensi standar yang punya mekanisme back-pressure.

Hasilnya cukup dramatis:

- Lead time rata-rata naik 4× (1.98 → 7.91)
- Persentil ke-95 lead time, worst-case yang dialami 5 persen order terparah, naik 6.6× (2.00 → 13.13)
- Variabilitas lead time naik 46×
- WIP total naik 3.7×
- Throughput turun 50%


╔══════════════════════════════════════════════════════════╗
║  UPLOAD GAMBAR  ►  results/07_bottleneck_effect.png       ║
║  Caption: Dampak bottleneck kapasitas: WIP jenuh di hulu, ║
║  distribusi lead time berekor panjang, umur WIP melonjak  ║
║  tepat di depan bottleneck.                               ║
╚══════════════════════════════════════════════════════════╝

Tiga hal yang terlihat jelas di gambar di atas. Pertama, WIP menumpuk jenuh di hulu bottleneck, buffer penuh karena barang tidak bisa lewat. Kedua, distribusi lead time berubah dari paku tajam yang bisa diprediksi jadi ekor panjang penuh ketidakpastian; manajer pabrik paling takut sama ekor panjang, karena artinya kadang-kadang order butuh waktu jauh di luar prediksi normal. Ketiga, umur WIP melonjak tepat di depan bottleneck, batch-batch di sana benar-benar terjebak.

Cerita ini kuantitatif tentang sesuatu yang orang produksi rasakan sehari-hari: satu stasiun yang lambat bisa melumpuhkan seluruh pabrik, dan gejala paling awalnya sering muncul jauh di hulu, bukan di titik masalahnya.

## Penutup, dan sedikit refleksi

Yang aku sukai dari proyek ini adalah betapa satu kerangka matematika, hukum kekekalan hiperbolik yang lahir dari dinamika gas dan lalu lintas, bisa bicara dengan fasih tentang lantai pabrik. Dan bagaimana dua sudut pandang yang tadinya kupikir bersaing, Eulerian dan Lagrangian, ternyata saling melengkapi kalau perannya diberi porsi yang tepat.

Pelajaran terbesar buatku bukan soal matematikanya, tapi soal kejujuran. Asumsi awalku soal Lagrangian keliru, dan mengakuinya justru membuka arah kontribusi yang lebih menarik dari rencana awal. Aku juga belajar bahwa mengerjakan riset dengan agen AI itu terasa kayak punya teman ngoprek yang selalu siap, selalu sabar, dan tidak segan kasih tahu kalau aku salah, bukan mengganti peran manusia, tapi memperbesar apa yang bisa dikerjakan satu manusia di waktu yang punya batas.

Proposal S2 yang lima tahun terkubur, akhirnya hidup lagi. Kali ini lebih sederhana, lebih dekat ke kerjaan sehari-hari, dan yang paling penting: ada kodenya yang beneran jalan.

Kode lengkap, semua eksperimen, dan proposalnya terbuka di: github.com/kandref/supply-chain-research

Kalau ada yang ingin ngobrol soal pemodelan numerik, rantai pasok, atau soal ngoprek bareng agen AI, silakan sapa.

---

_Ditulis oleh Kurnia Andre Febrian, alumnus Fisika Teoretis yang sekarang berlabuh di dunia BI. Solver dan plotnya dibangun dengan Python, NumPy, dan Matplotlib. Bagian iterasi cepatnya dibantu Hermes Agent (buatan Nous Research)._
