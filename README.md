# Komputasi Numerik: Penyelesaian Sistem Persamaan Linier (SPL)

Repositori ini berisi implementasi skrip komputasi numerik dalam bahasa MATLAB / GNU Octave untuk menyelesaikan Sistem Persamaan Linier (SPL) multivariabel (5x5). Proyek ini membandingkan efisiensi, akurasi, dan waktu komputasi dari lima algoritma penyelesaian matriks yang berbeda.

## 🚀 Fitur Utama

Proyek ini terbagi menjadi tiga blok komputasi utama:
1. **Penyelesaian Akar SPL**: Mencari nilai dari 5 variabel ($p, q, r, s, t$) menggunakan 5 metode numerik:
   - **Metode Langsung (Analitik):**
     * Eliminasi Gauss
     * Eliminasi Gauss-Jordan
     * Dekomposisi LU
   - **Metode Hampiran (Iteratif):**
     * Metode Jacobi
     * Metode Gauss-Seidel
2. **Analisis Galat (Error) Abstrak**: Menghitung selisih/deviasi nilai hampiran dari metode iteratif terhadap nilai eksak (Metode Gauss) dengan batasan toleransi error `0.001` dan maksimum iterasi `20`.
3. **Profiling Running Time**: Mengukur durasi waktu eksekusi murni dari tiap algoritma (algoritma kompleksitas) menggunakan fungsi `tic` dan `toc`, dengan mengeksklusi waktu latensi Input/Output (I/O).

## 🛠️ Persyaratan Sistem (Prerequisites)

Untuk menjalankan kode dalam repositori ini, Anda memerlukan salah satu dari perangkat lunak berikut:
- **MATLAB** (R2014a atau versi lebih baru)
- **GNU Octave** (Versi 5.x atau lebih baru) - *Gratis/Open Source*

## ⚙️ Cara Penggunaan (How to Use)

1. **Clone repositori ini** ke direktori lokal komputer Anda:
   ```bash
   git clone [https://github.com/username-anda/nama-repositori.git](https://github.com/username-anda/nama-repositori.git)
