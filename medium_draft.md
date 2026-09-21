# Ketika macet di jalan tol mengajari kita soal pabrik

*Cerita tentang memodelkan bottleneck rantai pasok dengan persamaan gelombang, dan tentang menghidupkan kembali sebuah proposal yang lima tahun terkubur.*

---

## Sedikit cerita, sebelum masuk ke matematika

Tahun 2021 aku sempat hampir mendaftar S2 Fisika Teoretis di ITERA. Judul proposalnya panjang dan agak sok penting, kira-kira "Analisis Penggunaan CFD untuk Pemetaan dan Rekayasa Lalu Lintas". Aku bikin sampai cover, lalu keburu dapat kerja dan proposalnya masuk laci. Selesai.

Lima tahun kemudian, di tengah kerjaan sehari-hari sebagai analis BI, aku iseng buka file lama itu. Tiba-tiba pengin coba lagi, tapi dengan sudut pandang yang lebih dekat ke kerjaan sekarang. Rantai pasok, bukan lalu lintas. Bottleneck produksi, bukan macet tol. Meskipun secara matematika, jujurnya, dua hal itu adalah adik-kakak.

Yang menarik, reset ini aku kerjakan bareng sebuah agen AI bernama Hermes Agent, sebuah agent framework buatan Nous Research yang bisa nulis kode, jalanin eksperimen numerik, compile LaTeX, dan pada saat yang sama juga jujur waktu asumsiku ternyata keliru (dan aku akan cerita soal itu di bawah). Bukan sekadar auto-complete. Lebih mirip teman ngoprek yang tahan begadang. Aku dikte arah dan justifikasi ilmiahnya, dia bantu turunkan jadi kode, plot, dan tulisan. Iterasi jadi cepat sekali. Yang biasanya butuh berminggu-minggu buat setup solver dan verifikasi, selesai dalam beberapa sesi obrolan.

Aku tulis ini sebagian sebagai catatan teknis, sebagian sebagai catatan pribadi. Kalau kamu tertarik ke sisi teknisnya saja, lompati bagian ini dan langsung ke bawah.

---

## Kenapa pabrik itu mirip jalan tol

Pernah lihat macet di tol yang bergerak mundur? Mobil di depan sudah jalan, tapi gelombang rem menyala terus merambat ke belakang, menabrak mobil yang baru datang. Fenomena ini punya nama, gelombang kejut lalu lintas, dan sudah lama dimodelkan dengan matematika yang sama yang dipakai buat aliran fluida.

Kejutannya, struktur matematis yang sama muncul juga di lantai pabrik. Ganti "mobil" dengan "produk setengah jadi", ganti "posisi di jalan" dengan "tahap penyelesaian produk", dan tiba-tiba kemacetan tol berubah jadi tumpukan barang di depan sebuah mesin yang kapasitasnya terbatas.

Alih-alih melacak jutaan unit produk satu per satu, kita perlakukan aliran produksi itu seperti fluida yang punya kerapatan ρ(x, t). Berapa banyak barang dalam proses yang ada pada tahap x pada waktu t. Sumbu x di sini bukan posisi fisik. Ia adalah "derajat penyelesaian" dari 0 (bahan mentah) sampai 1 (produk jadi).

Dinamikanya diatur satu persamaan sederhana:

> ∂ρ/∂t + ∂F(ρ)/∂x = 0

F(ρ) adalah *clearing function*, laju keluaran produksi sebagai fungsi seberapa padat WIP. Aku pakai bentuk yang populer di literatur rantai pasok, model Armbruster, Degond, Ringhofer:

> F(ρ) = μ · ρ / (1 + ρ)

Kalau WIP sedikit, output naik sebanding beban. Kalau WIP menumpuk, output jenuh ke kapasitas maksimum μ. Ini persis perilaku mesin nyata yang kadang bisa cepat, kadang macet.

---

## Membuktikan solvernya tidak bohong

Sebelum percaya hasil apa pun, solver numerik harus dites. Ini bagian yang sering dilewati orang, padahal paling penting. Aku sendiri nyaris melewatinya sebelum diingatkan.

*Riemann problem* adalah batu uji klasik. Dua kepadatan konstan bertemu di satu titik. Teori (kondisi Rankine dan Hugoniot) sudah kasih tahu persis seberapa cepat gelombang kejut harus merambat. Solverku (skema *finite volume* Godunov) menghasilkan posisi shock yang meleset 0.0006 dari nilai analitis. Kecil sekali. Boleh lega.

![Riemann problem: pembentukan gelombang kejut dan penjarangan](IMG_01_riemann)

Diagram ruang-waktu bikin fenomenanya lebih jelas. Gelombang kejut merambat sebagai garis diagonal yang tajam, persis analog dengan gelombang rem di tol tadi.

![Diagram ruang-waktu: bottleneck merambat dengan kecepatan Rankine-Hugoniot](IMG_01_spacetime)

Uji kedua, konvergensi. Kalau grid diperhalus, galat harus mengecil dengan laju yang bisa diprediksi. Hasilnya orde 1.0, persis sesuai teori skema Godunov orde satu. Sesuai teori dalam matematika numerik itu perasaan yang menyenangkan.

![Uji konvergensi: galat L1 mengecil sesuai orde teoretis](IMG_02_convergence)

Uji ketiga, konservasi massa. Barang tidak boleh hilang atau muncul dari ketiadaan. Residu neraca massaku 7 kali sepuluh pangkat minus delapan belas. Itu setara batas presisi mesin. Solver lolos.

---

## Dunia itu tidak pernah pasti

Pabrik nyata tidak deterministik. Order berfluktuasi, kapasitas mesin bergoyang, operator sakit, bahan baku telat. Kalau aku cuma jalankan simulasi sekali, hasilnya kelihatan meyakinkan padahal cuma satu skenario dari ribuan kemungkinan.

Jadi aku jalankan 400 kali dengan parameter acak. Kedatangan order pakai distribusi Poisson, kapasitas pakai normal. Hasilnya bukan satu angka throughput, tapi seluruh distribusi. Ini yang bikin model jadi berguna buat pengambil keputusan, karena mereka bisa lihat bukan cuma "rata-rata" tapi juga "seberapa buruk kemungkinan terburuknya".

![Distribusi throughput dari 400 simulasi Monte-Carlo](IMG_03_montecarlo)

---

## Dua sudut pandang, dan pengakuan bahwa aku salah

Ada dua cara memandang aliran. *Eulerian*, di mana kita duduk diam di satu titik dan mengamati kepadatan yang lewat. *Lagrangian*, di mana kita naik ke atas tiap paket produk dan ikut jalan bersamanya. Aku bangun keduanya dan cocokkan.

Pada kondisi mulus, keduanya sepakat sampai galat 2 persen. Bagus.

![Perbandingan silang Eulerian vs Lagrangian pada kondisi mulus](IMG_04_xcheck)

Lalu aku berasumsi, dengan cukup percaya diri, bahwa deskripsi Lagrangian akan lebih tajam di sekitar bottleneck. Alasannya intuitif: kalau kita ikut naik ke tiap paket, kita tahu persis di mana tiap paket berada, jadi transisi tajam mestinya lebih mudah tertangkap.

Ternyata aku salah.

![Coupling di shock: Godunov justru lebih tajam, Lagrangian overshoot](IMG_05_coupling)

Skema Godunov menghasilkan shock dengan lebar transisi 0.005, sementara rekonstruksi Lagrangian malah overshoot sampai kerapatan 3.8, dengan osilasi. Ini bukan bug yang bisa ditambal. Sifat teoretisnya memang begitu. Godunov punya properti monotonisitas yang bikin dia unggul di diskontinuitas, dan kernel rekonstruksi Lagrangian punya kelemahan alami di titik transisi tajam.

Alih-alih memaksakan asumsi lama, aku rumuskan ulang perannya. Lagrangian tidak untuk menajamkan shock. Lagrangian untuk melacak sesuatu yang Eulerian secara struktural tidak bisa: identitas dan riwayat tiap batch.

Momen ini yang paling berkesan buat aku selama ngerjain proyek ini. Waktu Hermes bantu diagnosa hasil yang aneh itu, dia tidak berusaha buat asumsi awalku terlihat benar. Dia bilang, dengan bahasa yang lebih halus, bahwa premis awalku keliru dan justru mengarahkan aku ke reframing yang jauh lebih menarik. Kejujuran seperti ini yang aku butuhkan.

---

## Umur dan lead time, di sinilah Lagrangian menemukan tempatnya

Deskripsi Eulerian tahu berapa banyak WIP di tiap tahap, tapi tidak tahu sudah berapa lama sebuah batch tertentu berada di sistem. Padahal itulah metrik yang manajer pabrik pedulikan setiap hari. Lead time. Umur produk. Umur bahan baku sebelum terpakai.

Dengan menandai tiap paket dengan waktu lahirnya, aku bisa merekam distribusi lead time produk yang selesai. Rata-ratanya keluar di 1.4914, sedangkan prediksi analitisnya 1.5000. Meleset 0.6 persen. Lumayan.

Tapi validasi yang paling meyakinkan datang dari hukum Little. Rumus ini bilang jumlah WIP di sistem sama dengan laju throughput dikali rata-rata lead time. Hukum ini berlaku universal untuk sistem antrian tunak apa pun, dari kasir minimarket sampai kilang minyak. Rasio prediksi hukum Little terhadap WIP aktual di simulasiku, 0.994. Meleset kurang dari 1 persen. Karena hukum ini begitu umum, kesesuaian ini adalah bukti kuat bahwa pelacakannya benar secara fisis.

![Pelacakan umur dan lead time: profil WIP, distribusi lead time, profil umur](IMG_06_aging)

---

## Menciptakan bottleneck dan menyaksikan kekacauan

Sekarang bagian paling seru. Aku turunkan kapasitas satu tahap ke 35 persen pada jendela sempit di x antara 0.55 sampai 0.65. Ini simulasi dari satu stasiun kerja yang kewalahan. Buat menghasilkan antrian yang merambat mundur (persis fenomena rem tol tadi), aku perluas flux ke bentuk *Cell Transmission Model*, ekstensi standar yang punya mekanisme *back-pressure*.

Hasilnya cukup dramatis. Lead time rata-rata naik empat kali lipat, dari 1.98 ke 7.91. Persentil ke-95 lead time (artinya, worst-case yang dialami 5 persen order terparah) naik 6.6 kali, dari 2.00 ke 13.13. Variabilitas lead time naik 46 kali. WIP total naik 3.7 kali. Throughput turun setengahnya.

![Dampak bottleneck: WIP menumpuk, lead time berekor panjang, umur melonjak](IMG_07_bottleneck)

Perhatikan tiga hal di gambar. Pertama, WIP menumpuk jenuh di hulu bottleneck. Buffer penuh karena barang tidak bisa lewat. Kedua, distribusi lead time berubah dari paku tajam yang bisa diprediksi jadi ekor panjang penuh ketidakpastian. Manajer pabrik paling takut sama ekor panjang, karena artinya kadang-kadang order butuh waktu jauh di luar prediksi normal. Ketiga, umur WIP melonjak tepat di depan bottleneck. Batch-batch di sana benar-benar terjebak.

Cerita ini kuantitatif tentang sesuatu yang orang produksi rasakan sehari-hari. Satu stasiun yang lambat bisa melumpuhkan seluruh pabrik, dan gejala paling awalnya sering muncul jauh di hulu, bukan di titik masalahnya.

---

## Penutup, dan sedikit refleksi

Yang aku sukai dari proyek ini adalah betapa satu kerangka matematika, hukum kekekalan hiperbolik yang lahir dari dinamika gas dan lalu lintas, bisa bicara dengan fasih tentang lantai pabrik. Dan bagaimana dua sudut pandang yang tadinya kupikir bersaing (Eulerian dan Lagrangian) ternyata saling melengkapi kalau perannya diberi porsi yang tepat.

Pelajaran terbesar buatku bukan soal matematikanya. Soal kejujuran. Asumsi awalku soal Lagrangian keliru, dan mengakuinya justru membuka arah kontribusi yang lebih menarik dari rencana awal. Aku juga belajar bahwa mengerjakan riset dengan agen AI itu terasa kayak punya teman ngoprek yang selalu siap, selalu sabar, dan tidak segan kasih tahu kalau aku salah. Bukan mengganti peran manusia, tapi memperbesar apa yang bisa dikerjakan satu manusia di waktu yang punya batas.

Proposal S2 yang lima tahun terkubur, akhirnya hidup lagi. Kali ini lebih sederhana, lebih dekat ke kerjaan sehari-hari, dan yang paling penting, ada kodenya yang beneran jalan.

Kode lengkap, semua eksperimen, dan proposalnya terbuka di:
**https://github.com/kandref/supply-chain-research**

Kalau ada yang ingin ngobrol soal pemodelan numerik, rantai pasok, atau soal ngoprek bareng agen AI, silakan sapa.

---

*Ditulis oleh Kurnia Andre Febrian, bekas anak Fisika Teoretis ITERA yang sekarang kerja di dunia BI. Solver dan plotnya dibangun dengan Python, NumPy, dan Matplotlib. Bagian iterasi cepatnya dibantu Hermes Agent (buatan Nous Research).*
