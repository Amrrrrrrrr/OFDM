close all; clc;

%% ========================================================================
%  TP TS205 - Partie 1 : OFDM sur canal AWGN complexe
%  Paramètres de simulation
%% ========================================================================
Nt  = 500;        % Nombre de symboles OFDM par trame
N   = 128;        % Nombre de sous-porteuses totales
Nu  = 128;        % Nombre de sous-porteuses utilisées
TCP = 0;          % Taille du préfixe cyclique (0 pour la partie 1)
Ts  = 0.05e-6;    % Temps symbole (0.05 µs)
K   = Nu * Nt;    % Nombre total de bits = 64000

%% ========================================================================
%  Question 1 : Validation de la chaîne (TEB = 0 sans bruit)
%% ========================================================================

% Source binaire aléatoire
b = randi([0 1], 1, K);

% Émetteur : BPSK + OFDM
x_mod = mod_BPSK(b);
signal_ofdm = ofdm_modulate(x_mod, N, Nu, Nt, TCP);

% Canal AWGN (sigma2 = 0 : pas de bruit)
signal_rx = canal_AWGN(signal_ofdm, 0);

% Récepteur : OFDM + BPSK
symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP);
b_r = demod_BPSK(symboles_rx);

% Vérification
nb_erreurs = sum(b ~= b_r);
TEB = nb_erreurs / K;

fprintf('============================================\n');
fprintf('  Question 1 : Validation de la chaîne\n');
fprintf('============================================\n');
fprintf('  Bits transmis : %d | Erreurs : %d | TEB : %e\n', K, nb_erreurs, TEB);
if TEB == 0
    fprintf('  >>> Test réussi ! TEB = 0 sans bruit.\n');
else
    fprintf('  >>> Test échoué ! Vérifier le code.\n');
end
fprintf('============================================\n\n');

% Constellation
figure;
scatter(real(symboles_rx), imag(symboles_rx), 1, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Partie réelle'); ylabel('Partie imaginaire');
title('Constellation BPSK après démodulation OFDM (sans bruit)');
grid on;

%% ========================================================================
%  Question 2 : Histogrammes des parties réelle et imaginaire
%% ========================================================================

signal_re = real(signal_ofdm);
signal_im = imag(signal_ofdm);

figure('Name', 'Histogrammes OFDM', 'Position', [100 100 1400 500]);

% Partie réelle
subplot(1, 2, 1);
histogram(signal_re, 100, 'FaceColor', 'b', 'Normalization', 'pdf');
hold on;
x_axis = linspace(min(signal_re), max(signal_re), 500);
mu_re = mean(signal_re);
sigma_re = std(signal_re);
gauss_re = (1/(sigma_re*sqrt(2*pi))) * exp(-0.5*((x_axis - mu_re)/sigma_re).^2);
plot(x_axis, gauss_re, 'r', 'LineWidth', 2);
hold off;
xlabel('Valeurs de la partie réelle', 'FontSize', 14);
ylabel('Densité de probabilité', 'FontSize', 14);
title('Partie réelle du signal OFDM', 'FontSize', 14);
legend('Histogramme', sprintf('Gaussienne N(%.2f, %.2f)', mu_re, sigma_re^2));
grid on;

% Partie imaginaire
subplot(1, 2, 2);
histogram(signal_im, 100, 'FaceColor', 'r', 'Normalization', 'pdf');
hold on;
mu_im = mean(signal_im);
sigma_im = std(signal_im);
x_axis = linspace(min(signal_im), max(signal_im), 500);
gauss_im = (1/(sigma_im*sqrt(2*pi))) * exp(-0.5*((x_axis - mu_im)/sigma_im).^2);
plot(x_axis, gauss_im, 'r', 'LineWidth', 2);
hold off;
xlabel('Valeurs de la partie imaginaire', 'FontSize', 14);
ylabel('Densité de probabilité', 'FontSize', 14);
title('Partie imaginaire du signal OFDM', 'FontSize', 14);
legend('Histogramme', sprintf('Gaussienne N(%.2f, %.2f)', mu_im, sigma_im^2));
grid on;

fprintf('============================================\n');
fprintf('  Question 2 : Statistiques du signal OFDM\n');
fprintf('============================================\n');
fprintf('  Re : moyenne = %.4f, variance = %.4f\n', mu_re, sigma_re^2);
fprintf('  Im : moyenne = %.4f, variance = %.4f\n', mu_im, sigma_im^2);
fprintf('  => Distribution gaussienne (TCL, N=%d)\n', N);
fprintf('============================================\n\n');

%% ========================================================================
%  Question 3 : Autocorrélation et intercorrélation
%% ========================================================================

autocorr_re = xcorr(signal_re, 'coeff');
autocorr_im = xcorr(signal_im, 'coeff');
intercorr   = xcorr(signal_re, signal_im, 'coeff');

L_sig  = length(signal_re);
lags   = -(L_sig-1):(L_sig-1);
center = L_sig;
idx = 1:length(lags);

figure('Name', 'Autocorrélation et intercorrélation', 'Position', [100 100 1600 450]);

subplot(1, 3, 1);
plot(lags(idx), autocorr_re(idx), 'b', 'LineWidth', 1.5);
xlabel('Lag (échantillons)', 'FontSize', 14);
ylabel('Autocorrélation', 'FontSize', 14);
title('Autocorrélation partie réelle', 'FontSize', 14);
grid on;

subplot(1, 3, 2);
plot(lags(idx), autocorr_im(idx), 'r', 'LineWidth', 1.5);
xlabel('Lag (échantillons)', 'FontSize', 14);
ylabel('Autocorrélation', 'FontSize', 14);
title('Autocorrélation partie imaginaire', 'FontSize', 14);
grid on;

subplot(1, 3, 3);
plot(lags(idx), intercorr(idx), 'g', 'LineWidth', 1.5);
xlabel('Lag (échantillons)', 'FontSize', 14);
ylabel('Intercorrélation', 'FontSize', 14);
title('Intercorrélation Re / Im', 'FontSize', 14);
grid on;

fprintf('============================================\n');
fprintf('  Question 3 : Corrélations\n');
fprintf('============================================\n');
fprintf('  - Autocorr : pic à lag=0, quasi nul ailleurs\n');
fprintf('    => signal quasi-blanc.\n');
fprintf('  - Périodicité tous les N=%d échantillons.\n', N);
fprintf('  - Intercorr Re/Im ≈ 0 => décorrélées.\n');
fprintf('============================================\n\n');

%% ========================================================================
%  Question 5 : DSP par périodogramme de Welch
%% ========================================================================

NFFT1 = N;
NFFT2 = 16 * N;
DSP_128  = Mon_Welch(signal_ofdm, NFFT1);
DSP_2048 = Mon_Welch(signal_ofdm, NFFT2);

figure('Name', 'DSP OFDM', 'Position', [100 100 1400 500]);

subplot(1, 2, 1);
freq1 = (0:NFFT1-1) / NFFT1;
plot(freq1, 10*log10(DSP_128 + 1e-12), 'b', 'LineWidth', 2);
xlabel('Fréquences normalisées', 'FontSize', 14);
ylabel('DSP (dB)', 'FontSize', 14);
title(sprintf('NFFT = %d', NFFT1), 'FontSize', 14);
grid on;

subplot(1, 2, 2);
freq2 = (0:NFFT2-1) / NFFT2;
plot(freq2, 10*log10(DSP_2048 + 1e-12), 'r', 'LineWidth', 1.5);
xlabel('Fréquences normalisées', 'FontSize', 14);
ylabel('DSP (dB)', 'FontSize', 14);
title(sprintf('NFFT = %d', NFFT2), 'FontSize', 14);
grid on;

fprintf('============================================\n');
fprintf('  Question 5 : DSP du signal OFDM\n');
fprintf('============================================\n');
fprintf('  NFFT = %d  : DSP plate (résolution = Df)\n', NFFT1);
fprintf('  NFFT = %d : sinc² visibles + lobes secondaires\n', NFFT2);
fprintf('============================================\n\n');

%% ========================================================================
%  Question 6 : TEB vs SNR (Monte Carlo)
%% ========================================================================

SNR_dB   = 0:1:10;
SNR_lin  = 10.^(SNR_dB / 10);
sigma2_S = 1;
var_bruit = sigma2_S ./ SNR_lin;
nb_errors_max = 100;

BER = zeros(1, length(SNR_dB));

for j = 1:length(SNR_dB)
    nb_bits_errors = 0;
    nb_bits_total  = 0;

    while nb_bits_errors < nb_errors_max * Nu
        b = randi([0 1], 1, K);
        x_mod     = mod_BPSK(b);
        signal_tx = ofdm_modulate(x_mod, N, Nu, Nt, TCP);
        signal_rx = canal_AWGN(signal_tx, var_bruit(j));
        symboles_rx = ofdm_demodulate(signal_rx, N, Nu, Nt, TCP);
        b_r = demod_BPSK(symboles_rx);

        nb_bits_errors = nb_bits_errors + sum(b ~= b_r);
        nb_bits_total  = nb_bits_total + K;
    end

    BER(j) = nb_bits_errors / nb_bits_total;
    fprintf('SNR = %2d dB : BER = %.2e (%d bits)\n', SNR_dB(j), BER(j), nb_bits_total);
end

% Théorie BPSK/AWGN : Pb = 0.5 * erfc(sqrt(SNR))
Pb_theo = 0.5 * erfc(sqrt(SNR_lin));

figure('Name', 'TEB OFDM vs BPSK théorique', 'Position', [100 100 700 500]);
semilogy(SNR_dB, BER,     'ro-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(SNR_dB, Pb_theo, 'b*-', 'LineWidth', 2, 'MarkerSize', 8); hold off;
xlabel('SNR (dB)', 'FontSize', 14);
ylabel('TEB', 'FontSize', 14);
title('TEB OFDM simulé vs BPSK théorique (canal AWGN)', 'FontSize', 14);
legend('OFDM simulé', 'BPSK théorique', 'FontSize', 12, 'Location', 'southwest');
grid on;

saveas(gcf, 'TEB_OFDM_vs_BPSK_AWGN.pdf');

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