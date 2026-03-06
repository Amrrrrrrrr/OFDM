close all; clc;

%% Paramètres
Nt  = 500;
N   = 128;
Nu  = 128;
TCP = 0;
Ts  = 0.05e-6;
K   = Nu * Nt;
%% ========================================================================
%  Partie 2 – Question 1 : Réponse fréquentielle du canal de Rayleigh
%% ========================================================================

list_L = [4, 8, 16, 32];
list_color = {'b', 'g', 'r', 'm'};
list_h_ray = cell(1, length(list_L));

figure('Name', 'Module du spectre du canal de Rayleigh');
hold on;

for l = 1:length(list_L)
    L = list_L(l);
    % Génération des coefficients : h ~ NC(0, 1/L)
    h = sqrt(1/(2*L)) * (randn(1, L) + 1j*randn(1, L));
    list_h_ray{l} = h;
    
    % FFT sur N points
    H = fft(h, N);
    
    plot(0:N-1, abs(H), list_color{l}, 'LineWidth', 1.5);
end

set(gca, 'YScale', 'log');
xlabel('Fréquences normalisées', 'FontSize', 14);
ylabel('|H[n]|', 'FontSize', 14);
title('Module du spectre du canal de Rayleigh', 'FontSize', 14);
legend('L = 4', 'L = 8', 'L = 16', 'L = 32', 'FontSize', 12);
grid on;
hold off;

%% ========================================================================
%  Partie 2 – Question 2 : Validation chaîne CP-OFDM
%% ========================================================================

TCP = 16;
K = Nu * Nt;

% Source binaire
b = randi([0 1], 1, K);

% Émetteur
x_mod = mod_BPSK(b);
signal_ofdm = ofdm_modulate(x_mod, N, Nu, Nt, TCP);

% Canal trivial : h = delta[p], pas de bruit
signal_rx = canal_AWGN(signal_ofdm, 0);

% Récepteur
symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP);
b_r = demod_BPSK(symboles_rx);

% Vérification
nb_erreurs = sum(b ~= b_r);
if nb_erreurs == 0
    fprintf('Partie 2 - Q2 : Erreurs = %d => Test réussi!\n', nb_erreurs);
else
    fprintf('Partie 2 - Q2 : Erreurs = %d => Test échoué!\n', nb_erreurs);
end

%% ========================================================================
%  Partie 2 – Question 4 : Égaliseur Zero Forcing
%% ========================================================================

TCP = 16;
K = Nu * Nt;

% Canal de Rayleigh avec L = 16
h = list_h_ray{3};  % L = 16 (3ème élément)

% Réponse fréquentielle du canal et son inverse
H = fft(h, N);
H_inv = 1 ./ H;

% Source binaire
b = randi([0 1], 1, K);

% Émetteur
x_mod = mod_BPSK(b);
signal_ofdm = ofdm_modulate(x_mod, N, Nu, Nt, TCP);

% Canal de Rayleigh (sans bruit)
signal_rx = canal_Rayleigh(signal_ofdm, h);

% Démodulation OFDM
symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP);

% Égalisation ZF
symboles_eq = EQ_ZF(symboles_rx, H_inv, Nu, Nt);

% Démodulation BPSK
b_r = demod_BPSK(symboles_eq);

% Vérification
nb_erreurs = sum(b ~= b_r);
if nb_erreurs == 0
    fprintf('Partie 2 - Q4 : Erreurs = %d => Test réussi!\n', nb_erreurs);
else
    fprintf('Partie 2 - Q4 : Erreurs = %d => Test échoué!\n', nb_erreurs);
end

% Constellation avant/après égalisation (sous-porteuse 42)
n_sp = 42;
symb_mat_rx = reshape(symboles_rx, Nu, Nt);
symb_mat_eq = reshape(symboles_eq, Nu, Nt);

figure;
scatter(real(symb_mat_rx(n_sp,:)), imag(symb_mat_rx(n_sp,:)), 10, 'r', 'filled'); hold on;
scatter(real(symb_mat_eq(n_sp,:)), imag(symb_mat_eq(n_sp,:)), 10, 'b', 'filled');
xlabel('Partie réelle', 'FontSize', 14);
ylabel('Partie imaginaire', 'FontSize', 14);
title('Constellation sous-porteuse 42 : avant/après égalisation', 'FontSize', 14);
legend('Non égalisé', 'Égalisé', 'FontSize', 12);
grid on;

%% ========================================================================
%  Partie 2 – Question 5 : Module et Phase de H[n] pour L = 16
%% ========================================================================

h = list_h_ray{3};  % L = 16
H = fft(h, N);

figure('Name', 'Module et Phase du canal Rayleigh', 'Position', [100 100 1400 500]);

% Module
subplot(1, 2, 1);
plot(0:N-1, abs(H), 'b', 'LineWidth', 2);
xlabel('Fréquences normalisées', 'FontSize', 14);
ylabel('|H[n]|', 'FontSize', 14);
title('Module de H[n] (L=16)', 'FontSize', 14);
set(gca, 'YScale', 'log');
grid on;

% Phase
subplot(1, 2, 2);
plot(0:N-1, angle(H), 'r', 'LineWidth', 2);
xlabel('Fréquences normalisées', 'FontSize', 14);
ylabel('Phase (rad)', 'FontSize', 14);
title('Phase de H[n] (L=16)', 'FontSize', 14);
grid on;

%% ========================================================================
%  Partie 2 – Question 8 : Bandes de garde (30 sous-porteuses de chaque côté)
%% ========================================================================

Nu_guard = 128 - 2*30;  % 68 sous-porteuses utiles
TCP = 16;
K_guard = Nu_guard * Nt;

% Source
b = randi([0 1], 1, K_guard);
x_mod = mod_BPSK(b);

% Modulation OFDM avec bandes de garde
signal_ofdm_guard = ofdm_modulate_guard(x_mod, N, Nu_guard, Nt, TCP);

% DSP avec NFFT = 4*N
NFFT = 4 * N;
DSP_guard = Mon_Welch(signal_ofdm_guard, NFFT);

% Plot
figure('Name', 'DSP OFDM avec bandes de garde');
freq = (0:NFFT-1) / NFFT;
plot(freq, 10*log10(DSP_guard + 1e-12), 'b', 'LineWidth', 2);
xlabel('Fréquences normalisées', 'FontSize', 14);
ylabel('DSP (dB)', 'FontSize', 14);
title('DSP OFDM avec bandes de garde (30 de chaque côté)', 'FontSize', 14);
grid on;

%% ========================================================================
%  Partie 2 – Question 6 : TEB avec égalisation ZF sur canal Rayleigh
%% ========================================================================

h = list_h_ray{3};  % L = 16
H = fft(h, N);
H_inv = 1 ./ H;

TCP = 16;
K = Nu * Nt;

SNR_dB = 0:1:15;
SNR_lin = 10.^(SNR_dB / 10);
var_bruit = 1 ./ SNR_lin;   % sigma2_S = 1, E(|H|^2) = 1

nb_errors_max = 100;
BER_eq = zeros(1, length(SNR_dB));

for j = 1:length(SNR_dB)
    nb_bits_errors = 0;
    nb_bits_total = 0;

    while nb_bits_errors < nb_errors_max * Nu
        b = randi([0 1], 1, K);
        x_mod = mod_BPSK(b);
        signal_ofdm = ofdm_modulate(x_mod, N, Nu, Nt, TCP);

        % Canal Rayleigh + AWGN
        signal_ray = canal_Rayleigh(signal_ofdm, h);
        signal_rx = canal_AWGN(signal_ray, var_bruit(j));

        % Démodulation + Égalisation ZF
        symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP);
        symboles_eq = EQ_ZF(symboles_rx, H_inv, Nu, Nt);
        b_r = demod_BPSK(symboles_eq);

        nb_bits_errors = nb_bits_errors + sum(b ~= b_r);
        nb_bits_total = nb_bits_total + K;
    end

    BER_eq(j) = nb_bits_errors / nb_bits_total;
    fprintf('SNR = %2d dB : BER = %.2e\n', SNR_dB(j), BER_eq(j));
end

% Courbe théorique BPSK/AWGN
Pb_theo = 0.5 * erfc(sqrt(SNR_lin));

% Plot
figure('Name', 'TEB OFDM Rayleigh avec EQ ZF');
semilogy(SNR_dB, BER_eq, 'ro-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_dB, Pb_theo, 'b*-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)', 'FontSize', 14);
ylabel('TEB', 'FontSize', 14);
title('TEB OFDM + EQ ZF sur canal Rayleigh vs BPSK théorique', 'FontSize', 14);
legend('OFDM + ZF (Rayleigh)', 'BPSK théorique (AWGN)', 'FontSize', 12, 'Location', 'southwest');
grid on;

%% ========================================================================
%  Partie 2 – Question 7 : CP < L (CP = 8, L = 16)
%% ========================================================================

TCP_short = 8;
K = Nu * Nt;

h = list_h_ray{3};  % L = 16
H = fft(h, N);
H_inv = 1 ./ H;

SNR_dB = 0:1:15;
SNR_lin = 10.^(SNR_dB / 10);
var_bruit = 1 ./ SNR_lin;

nb_errors_max = 100;
BER_cp8 = zeros(1, length(SNR_dB));

for j = 1:length(SNR_dB)
    nb_bits_errors = 0;
    nb_bits_total = 0;

    while nb_bits_errors < nb_errors_max * Nu
        b = randi([0 1], 1, K);
        x_mod = mod_BPSK(b);
        signal_ofdm = ofdm_modulate(x_mod, N, Nu, Nt, TCP_short);

        signal_ray = canal_Rayleigh(signal_ofdm, h);
        signal_rx = canal_AWGN(signal_ray, var_bruit(j));

        symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP_short);
        symboles_eq = EQ_ZF(symboles_rx, H_inv, Nu, Nt);
        b_r = demod_BPSK(symboles_eq);

        nb_bits_errors = nb_bits_errors + sum(b ~= b_r);
        nb_bits_total = nb_bits_total + K;
    end

    BER_cp8(j) = nb_bits_errors / nb_bits_total;
    fprintf('SNR = %2d dB : BER (CP=8) = %.2e\n', SNR_dB(j), BER_cp8(j));
end

% Plot comparatif
Pb_theo = 0.5 * erfc(sqrt(SNR_lin));

figure('Name', 'TEB : CP=16 vs CP=8');
semilogy(SNR_dB, BER_eq,  'ro-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_dB, BER_cp8, 'gs-', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(SNR_dB, Pb_theo, 'b*-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('SNR (dB)', 'FontSize', 14);
ylabel('TEB', 'FontSize', 14);
title('Impact du CP insuffisant (L=16)', 'FontSize', 14);
legend('CP=16 (ZF)', 'CP=8 < L (ZF)', 'BPSK théorique (AWGN)', 'FontSize', 12, 'Location', 'southwest');
grid on;
%% ========================================================================
%  FONCTIONS
%% ========================================================================

function s = mod_BPSK(bits)
    s = 2 * bits - 1;
end

function bits = demod_BPSK(symboles)
    bits = double(real(symboles) > 0);
end

function signal_rx = canal_AWGN(signal_tx, sigma2)
    bruit = sqrt(sigma2/2) * (randn(size(signal_tx)) + 1j*randn(size(signal_tx)));
    signal_rx = signal_tx + bruit;
end

function signal_ofdm = ofdm_modulate(x_mod, N, Nu, Nt, TCP)
    X = reshape(x_mod, Nu, Nt);
    s = ifft(X, N, 1) * sqrt(N);
    if TCP > 0
        cp = s(end-TCP+1:end, :);
        s = [cp; s];
    end
    signal_ofdm = s(:).';
end

function symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP)
    r = reshape(signal_rx, N + TCP, Nt);
    if TCP > 0
        r = r(TCP+1:end, :);
    end
    R = fft(r, N, 1) / sqrt(N);
    R_util = R(1:Nu, :);
    symboles_rx = R_util(:).';
end

function DSP = Mon_Welch(signal, NFFT)
    L = length(signal);
    nb_segments = floor(L / NFFT);
    DSP = zeros(1, NFFT);
    for i = 0:nb_segments-1
        segment = signal(i*NFFT+1 : (i+1)*NFFT);
        S = fft(segment, NFFT) / NFFT;
        DSP = DSP + abs(S).^2;
    end
    DSP = DSP / nb_segments;
end

function y = canal_Rayleigh(signal, h)
    y = filter(h, 1, signal);
end

function symboles_eq = EQ_ZF(symboles_rx, H_inv, Nu, Nt)
    X = reshape(symboles_rx, Nu, Nt);
    H_inv_col = H_inv(1:Nu).';  % vecteur colonne
    X_eq = X .* H_inv_col;       % multiplication élément par élément
    symboles_eq = X_eq(:).';
end

function signal_ofdm = ofdm_modulate_guard(x_mod, N, Nu, Nt, TCP)
    % Nu = nombre de sous-porteuses utiles (68)
    % On place les données au centre, avec 30 zéros à gauche et à droite
    X_data = reshape(x_mod, Nu, Nt);
    X = zeros(N, Nt);
    % Sous-porteuses 31 à 98 portent les données (indices MATLAB : 31:98)
    X(31:30+Nu, :) = X_data;
    
    s = ifft(X, N, 1) * sqrt(N);
    if TCP > 0
        cp = s(end-TCP+1:end, :);
        s = [cp; s];
    end
    signal_ofdm = s(:).';
end