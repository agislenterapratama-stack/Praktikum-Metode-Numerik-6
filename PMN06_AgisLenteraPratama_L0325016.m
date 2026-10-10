% =========================================================================
% TUGAS PENYELESAIAN SISTEM PERSAMAAN LINIER
% Programmer : Agis Lentera Pratama
% NIM        : L0325016
% =========================================================================
clear; clc;

% --- INISIALISASI INPUT GLOBAL ---
% Matriks koefisien A dan vektor konstanta b dari persamaan tugasSPL.pdf[cite: 1]
A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];
b = [-1; 2; 6; 2; -1];

n_vars = length(b);
max_iter = 20;      % Maksimum iterasi iteratif
tol = 0.001;        % Batas toleransi error

22
%% ========================================================================
% SOAL 1: PENYELESAIAN SPL DENGAN 5 METODE
% =========================================================================
fprintf('================ SOAL 1: HASIL PENYELESAIAN ================\n');

% 1. Metode Gauss
Ag = [A b];
for i = 1:n_vars-1
    for j = i+1:n_vars
        m = Ag(j,i)/Ag(i,i);
        Ag(j,:) = Ag(j,:) - m*Ag(i,:);
    end
end
x_gauss = zeros(n_vars,1);
x_gauss(n_vars) = Ag(n_vars,n_vars+1)/Ag(n_vars,n_vars);
for i = n_vars-1:-1:1
    x_gauss(i) = (Ag(i,n_vars+1) - Ag(i,i+1:n_vars)*x_gauss(i+1:n_vars))/Ag(i,i);
end
disp(['1. Gauss        : ', num2str(x_gauss')]);

% 2. Metode Gauss-Jordan
Agj = [A b];
for i = 1:n_vars
    Agj(i,:) = Agj(i,:) / Agj(i,i);
    for j = 1:n_vars
        if i ~= j
            Agj(j,:) = Agj(j,:) - Agj(j,i)*Agj(i,:);
        end
    end
end
x_gj = Agj(:,n_vars+1);
disp(['2. Gauss-Jordan : ', num2str(x_gj')]);

% 3. Metode Dekomposisi LU
[L, U] = lu(A);
y = L \ b;
x_lu = U \ y;
disp(['3. LU Decomp    : ', num2str(x_lu')]);

% 4. Metode Jacobi
x_jac = zeros(n_vars,1);
x_new = zeros(n_vars,1);
for iter = 1:max_iter
    for i = 1:n_vars
        sigma = A(i,:)*x_jac - A(i,i)*x_jac(i);
        x_new(i) = (b(i) - sigma)/A(i,i);
    end
    if max(abs(x_new - x_jac)) < tol
        break;
    end
    x_jac = x_new;
end
disp(['4. Jacobi       : ', num2str(x_jac')]);

% 5. Metode Gauss-Seidel
x_gs = zeros(n_vars,1);
for iter = 1:max_iter
    x_old = x_gs;
    for i = 1:n_vars
        sigma = A(i,:)*x_gs - A(i,i)*x_gs(i);
        x_gs(i) = (b(i) - sigma)/A(i,i);
    end
    if max(abs(x_gs - x_old)) < tol
        break;
    end
end
disp(['5. Gauss-Seidel : ', num2str(x_gs')]);
fprintf('\n');


%% ========================================================================
% SOAL 2: PERHITUNGAN GALAT / ERROR
% =========================================================================
fprintf('================ SOAL 2: GALAT ABSTRAK =====================\n');
% Nilai eksak diambil dari metode analitik (Gauss)[cite: 2]
% Galat abstrak dihitung dengan rumus: abs(nilai eksak - nilai hampiran)[cite: 2]

galat_jacobi = abs(x_gauss - x_jac);
galat_gs = abs(x_gauss - x_gs);

disp(['Galat Metode Jacobi       : ', num2str(galat_jacobi')]);
disp(['Galat Metode Gauss-Seidel : ', num2str(galat_gs')]);
fprintf('\n');


%% ========================================================================
% SOAL 3: RUNNING TIME (WAKTU KOMPUTASI)
% =========================================================================
fprintf('================ SOAL 3: RUNNING TIME (DETIK) ==============\n');
% Pengukuran menggunakan tic toc murni pada blok algoritma tanpa I/O[cite: 2]

% 1. Waktu Gauss
tic;
Ag_rt = [A b];
for i = 1:n_vars-1
    for j = i+1:n_vars
        m = Ag_rt(j,i)/Ag_rt(i,i);
        Ag_rt(j,:) = Ag_rt(j,:) - m*Ag_rt(i,:);
    end
end
x_g_rt = zeros(n_vars,1);
x_g_rt(n_vars) = Ag_rt(n_vars,n_vars+1)/Ag_rt(n_vars,n_vars);
for i = n_vars-1:-1:1
    x_g_rt(i) = (Ag_rt(i,n_vars+1) - Ag_rt(i,i+1:n_vars)*x_g_rt(i+1:n_vars))/Ag_rt(i,i);
end
waktu_gauss = toc;
disp(['Waktu Gauss        : ', num2str(waktu_gauss)]);

% 2. Waktu Gauss-Jordan
tic;
Agj_rt = [A b];
for i = 1:n_vars
    Agj_rt(i,:) = Agj_rt(i,:) / Agj_rt(i,i);
    for j = 1:n_vars
        if i ~= j
            Agj_rt(j,:) = Agj_rt(j,:) - Agj_rt(j,i)*Agj_rt(i,:);
        end
    end
end
x_gj_rt = Agj_rt(:,n_vars+1);
waktu_gj = toc;
disp(['Waktu Gauss-Jordan : ', num2str(waktu_gj)]);

% 3. Waktu LU Decomp
tic;
[L_rt, U_rt] = lu(A);
y_rt = L_rt \ b;
x_lu_rt = U_rt \ y_rt;
waktu_lu = toc;
disp(['Waktu LU Decomp    : ', num2str(waktu_lu)]);

% 4. Waktu Jacobi
tic;
x_j_rt = zeros(n_vars,1);
x_n_rt = zeros(n_vars,1);
for iter = 1:max_iter
    for i = 1:n_vars
        sigma = A(i,:)*x_j_rt - A(i,i)*x_j_rt(i);
        x_n_rt(i) = (b(i) - sigma)/A(i,i);
    end
    if max(abs(x_n_rt - x_j_rt)) < tol
        break;
    end
    x_j_rt = x_n_rt;
end
waktu_jacobi = toc;
disp(['Waktu Jacobi       : ', num2str(waktu_jacobi)]);

% 5. Waktu Gauss-Seidel
tic;
x_gs_rt = zeros(n_vars,1);
for iter = 1:max_iter
    x_o_rt = x_gs_rt;
    for i = 1:n_vars
        sigma = A(i,:)*x_gs_rt - A(i,i)*x_gs_rt(i);
        x_gs_rt(i) = (b(i) - sigma)/A(i,i);
    end
    if max(abs(x_gs_rt - x_o_rt)) < tol
        break;
    end
end
waktu_gs = toc;
disp(['Waktu Gauss-Seidel : ', num2str(waktu_gs)]);
fprintf('============================================================\n');
