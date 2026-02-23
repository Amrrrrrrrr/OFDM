%% ========================================================================
%  TP TS205 - Partie 1 : OFDM sur canal AWGN complexe
%  Question 1 : Implémentation de la chaîne de communication OFDM
%  Validation : TEB = 0 lorsque sigma_nl^2 = 0
%% ========================================================================
clear all; close all; clc;

%% --- Paramètres de simulation -------------------------------------------
% Modulation
mod_order = 2;              % BPSK (M = 2)

% Temps symbole
Ts = 0.05e-6;               % Ts = 0.05 µs

% Nombre de sous-porteuses
N  = 128;                   % Nombre de sous-porteuses totales
Nu = 128;                   % Nombre de sous-porteuses utilisées

% Trame OFDM
NT = 500;                   % Nombre de symboles OFDM par trame

% Canal de propagation : h_l[p] = delta[p] (canal idéal)
h = 1;                      % Réponse impulsionnelle = impulsion de Dirac

% Variance du bruit (= 0 pour la validation)
sigma2_nl = 0;

%% --- Émetteur OFDM ------------------------------------------------------

% 1) Génération des bits aléatoires
%    Chaque symbole OFDM transporte Nu bits (BPSK : 1 bit/symbole)
nb_bits_total = Nu * NT;
bits_tx = randi([0 1], nb_bits_total, 1);

% 2) Modulation BPSK : 0 -> -1, 1 -> +1
symboles_tx = 2 * bits_tx - 1;

% 3) Conversion série -> parallèle : matrice Nu x NT
%    Chaque colonne = 1 symbole OFDM (N symboles sur N sous-porteuses)
S_tx = reshape(symboles_tx, Nu, NT);

% 4) Modulation OFDM : IFFT colonne par colonne
%    s_k = F^H * S_k  (IFFT normalisée par sqrt(N) pour conserver l'énergie)
s_tx = ifft(S_tx, N, 1) * sqrt(N);

%% --- Canal de propagation -----------------------------------------------

% 5) Convolution avec le canal h_l[p] = delta[p]
%    Ici h = 1, donc r_l[p] = s_l[p] (pas de distorsion)
s_canal = s_tx;             % Pas de filtrage (canal idéal)

% 6) Ajout du bruit AWGN complexe
%    n_l[p] ~ NC(0, sigma2_nl)
bruit = sqrt(sigma2_nl / 2) * (randn(N, NT) + 1j * randn(N, NT));
r_rx = s_canal + bruit;

%% --- Récepteur OFDM -----------------------------------------------------

% 7) Démodulation OFDM : FFT colonne par colonne
R_rx = fft(r_rx, N, 1) / sqrt(N);

% 8) Extraction des sous-porteuses utilisées
R_rx_util = R_rx(1:Nu, :);

% 9) Décision BPSK : seuil à 0 sur la partie réelle
symboles_rx = sign(real(R_rx_util));

% 10) Démodulation BPSK : -1 -> 0, +1 -> 1
bits_rx = (symboles_rx + 1) / 2;

%% --- Calcul du TEB ------------------------------------------------------

% Conversion en vecteur colonne
bits_rx_vect = reshape(bits_rx, nb_bits_total, 1);

% Calcul du taux d'erreur binaire
nb_erreurs = sum(bits_tx ~= bits_rx_vect);
TEB = nb_erreurs / nb_bits_total;

%% --- Affichage des résultats --------------------------------------------
fprintf('============================================\n');
fprintf('  TP TS205 - OFDM sur canal AWGN\n');
fprintf('  Question 1 : Validation de la chaîne\n');
fprintf('============================================\n');
fprintf('  Paramètres :\n');
fprintf('    Modulation       : BPSK\n');
fprintf('    N (sous-port.)   : %d\n', N);
fprintf('    Nu (utilisées)   : %d\n', Nu);
fprintf('    NT (symb. OFDM)  : %d\n', NT);
fprintf('    Ts               : %.2f µs\n', Ts * 1e6);
fprintf('    sigma2_nl        : %g\n', sigma2_nl);
fprintf('--------------------------------------------\n');
fprintf('  Résultats :\n');
fprintf('    Bits transmis    : %d\n', nb_bits_total);
fprintf('    Erreurs          : %d\n', nb_erreurs);
fprintf('    TEB              : %e\n', TEB);
fprintf('============================================\n');

if TEB == 0
    fprintf('  >>> VALIDATION OK : TEB = 0 (sans bruit)\n');
else
    fprintf('  >>> ERREUR : TEB non nul ! Vérifier le code.\n');
end
fprintf('============================================\n');