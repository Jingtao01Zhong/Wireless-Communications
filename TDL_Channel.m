%% Task 1 (MIMO: MF, ZF and MMSE)

%%
clear all;

%% Initialization (Feel free to experiment with these values)

% Number of Tx (transmitter) antennas
N_Tx = 2;

% Number of Rx (receiver) antennas
N_Rx = 4;

% SNR to evaluate
SNR_dB = 10;
SNR=10^(SNR_dB/10);

% Number of symbols per transmitter 
N_symbols = 1000;


%% Gray constellation for 16-QAM 
M = 16;
Constellation_16QAM = qammod(0:M-1, M, UnitAveragePower=true);

% Constellation scatter plot
figure(1); clf;
hold on
scatter(real(Constellation_16QAM(:)), imag(Constellation_16QAM(:)), 30,'filled')
line([-1.5 1.5], [0 0], 'Color', 'black')
line([0 0], [-1.5 1.5], 'Color', 'black')
xlabel("real"); ylabel("imag");
title('16 QAM constellation')
grid on; grid minor; axis square;


%% Data Deneration

% Generating random 16-QAM symbols
Orig_symbols = randi([0, M-1], N_Tx, N_symbols);
% Builtin QAM modulation function (qammod)
Mapped_symbols = qammod(Orig_symbols, M, UnitAveragePower=true);


%% Channel and Noise Realization

% The channel matrix H is assumed to be constant for all symbols.

% The (i,j) element of H is the channel between the 
% i:th RECEIVER antenna and the j:th TRANSMITTER antenna (y = Hx + n).

% Each of these channels is assumed to be an independent Rayleigh channel.
H = 1/sqrt(2) .* (randn(N_Rx, N_Tx) + 1i*randn(N_Rx, N_Tx));

% Average energy of received signal (signal power is 1: ||x|| = 1)
Es = mean(diag(H*H'));

% Generating noise to match the SNR 
sigma_Noise = sqrt(Es/SNR);
Noise = sigma_Noise * 1/sqrt(2) .* (randn(N_Rx, N_symbols) + 1i*randn(N_Rx, N_symbols));


%% Transmitting Data Through the Channel

% y = Hx + n for each group of N_tx transmitted symbols
Received_symbols = zeros(N_Rx, N_symbols);
for k = 1:N_symbols
    Received_symbols(:, k) = H * Mapped_symbols(:, k) + Noise(:, k);
end


%% Linear Receiver Implementations (MF, ZF and MMSE)

% Here we want to implement linear estimator matrices F such that
% x_est = Fy = FHx + Fn 
% can be used to estimate the sent symbol vector x

% Matched filter (MF) 
% (calling diag twice returns a matrix with only the diagonal elements left)
% This is basically just H' with some scaling for the rows
F_MF = inv(diag(diag(H'*H))) * H';


% Zero forcing (ZF) 
% F is set to be the inverse (or Moore Penrose pseudo-inverse) of H
% (inverting the channel)

% TASK 1A: FILL IN THE DIFFERENT CASES FOR ZERO FORCING FILTER

% Based on the number of antennas on Tx and Rx, either 
% H'*H (channel covariance) or H*H'(Rx channel covariance) is invertible
if (N_Rx == N_Tx)
    % Rectangular H, can be inverted directly
    F_ZF = H^-1;

elseif (N_Rx > N_Tx)
    % More receiver antennas (use pseudoinverse)
    F_ZF = (H' * H)^-1 * H';

else
    % More transmitter antennas (use pseudoinverse)
    F_ZF = H' * (H * H')^-1;

end


% Minimum mean square estimator (MMSE)
N_0 = sigma_Noise^2; % Noise power

% TASK 1B: FILL IN THE FORMULA FOR THE MMSE FILTER
% You can use eye(n) to generate a n X n identity matrix (n = N_Tx or N_Rx)
F_MMSE = (H' * H + N_0*eye(N_Tx))^-1 * H';
 

%% Filtering the Received Symbols

% Keeping track of the filtered symbols
Z_MF = zeros(N_Tx, N_symbols);
Z_ZF = zeros(N_Tx, N_symbols);
Z_MMSE = zeros(N_Tx, N_symbols);

% Loop through each group of received symbols y (= Hx + n)
for k = 1:N_symbols
    % The N_Tx x 1 vector of received symbols
    y = Received_symbols(:, k);
    
    % TASK 1C: USE THE FILTERS TO ESTIMATE THE ORIGINAL SIGNAL VECTOR x 
    % FROM THE RECEIVED SIGNAL VECTOR y

    % Matched filter (MF) (F_MF)
    Z_MF(:,k) = F_MF * y;

    % Zero forcing (ZF) (F_ZF)
    Z_ZF(:,k) = F_ZF * y;

    % Minimum mean square estimator (MMSE) (F_MMSE)
    Z_MMSE(:,k) = F_MMSE * y;
end

%% Scatter Plots for the Received and Filtered Symbols

% Received symbols
figure(2); clf;
subplot(2,2,1);
scatter(real(Received_symbols(:)), imag(Received_symbols(:)), 3,'filled')
title('Constellation of Received Symbols')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

% MF filter
subplot(2,2,2);
scatter(real(Z_MF(:)), imag(Z_MF(:)), 3,'filled')
title('Constellation of MF Receiver')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

% ZF filter
subplot(2,2,3);
scatter(real(Z_ZF(:)), imag(Z_ZF(:)), 3,'filled')
title('Constellation of ZF Receiver')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

% MMSE filter
subplot(2,2,4);
scatter(real(Z_MMSE(:)), imag(Z_MMSE(:)), 3,'filled')
title('Constellation of MMSE Receiver')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;


%% SINR of each Estimated Sample (Lec 6 Appendix)

% Signal and noise powers after filtering
Signal_power_MF = abs(F_MF*H).^2;
Noise_power_MF = N_0*diag(F_MF*F_MF');

Signal_power_ZF = abs(F_ZF*H).^2;
Noise_power_ZF = N_0*diag(F_ZF*F_ZF');

Signal_power_MMSE = abs(F_MMSE*H).^2;
Noise_power_MMSE = N_0*diag(F_MMSE*F_MMSE');

% Loop through the different transmitter antennas
SINR_vec_MF = zeros(1,N_Tx);
SINR_vec_ZF = zeros(1,N_Tx);
SINR_vec_MMSE = zeros(1,N_Tx);
for k = 1:N_Tx
    % Signal power
    S_k_MF = Signal_power_MF(k,k);
    S_k_ZF = Signal_power_ZF(k,k);
    S_k_MMSE = Signal_power_MMSE(k,k);

    % Noise+Interference power
    I_k_MF = sum(Signal_power_MF(k, :))-S_k_MF + Noise_power_MF(k);
    I_k_ZF = sum(Signal_power_ZF(k, :))-S_k_ZF + Noise_power_ZF(k);
    I_k_MMSE = sum(Signal_power_MMSE(k, :))-S_k_MMSE + Noise_power_MMSE(k);
   
    % SINR
    SINR_vec_MF(k) = S_k_MF/I_k_MF;
    SINR_vec_ZF(k) = S_k_ZF/I_k_ZF;
    SINR_vec_MMSE(k) = S_k_MMSE/I_k_MMSE;
end


% SINR plots
figure(3);clf;
hold on
plot(1:N_Tx, 10*log10(SINR_vec_MF))
plot(1:N_Tx, 10*log10(SINR_vec_ZF), '--', 'Linewidth', 2)
plot(1:N_Tx, 10*log10(SINR_vec_MMSE))
legend(["MF", "ZF", "MMSE"], 'Location', 'best')
title("SINR of the estimated signal values")
ylabel("SINR (dB)")
xlabel("Index of the transmitter antenna")
grid on; grid minor;


%% Symbol Error Rates

% Builtin QAM demodulation function (qamdemod)

% MF
Symbol_estimate_MF = qamdemod(Z_MF, M, UnitAveragePower=true);
Symbol_error_rate_MF = 1 - mean(Symbol_estimate_MF == Orig_symbols, 'all')

% ZF
Symbol_estimate_ZF = qamdemod(Z_ZF, M, UnitAveragePower=true);
Symbol_error_rate_ZF = 1 - mean(Symbol_estimate_ZF == Orig_symbols, 'all')

% MMSE
Symbol_estimate_MMSE = qamdemod(Z_MMSE, M, UnitAveragePower=true);
Symbol_error_rate_MMSE = 1 - mean(Symbol_estimate_MMSE == Orig_symbols, 'all')





