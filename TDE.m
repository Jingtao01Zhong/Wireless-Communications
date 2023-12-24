%% Task 2 (Time-Domain Equalization (TDE))

%%
clear all;

%% Initialization (Feel free to experiment with these values)


% Receiver block size that is used to estimate the transmitted signal
N_r = 5;

% SNR to evaluate
SNR_dB = 30;
SNR=10^(SNR_dB/10);


% Number of symbols
N_symbols = 1e5;

%% Gray constellation for 16-QAM 
M = 16;
Constellation_16QAM = qammod(0:M-1, M, UnitAveragePower=true);

% Constellation scatter plot
figure(4); clf;
hold on
scatter(real(Constellation_16QAM(:)), imag(Constellation_16QAM(:)), 30,'filled')
line([-1.5 1.5], [0 0], 'Color', 'black')
line([0 0], [-1.5 1.5], 'Color', 'black')
xlabel("real"); ylabel("imag");
title('16 QAM constellation')
grid on; grid minor; axis square;

%% Data Deneration

% Generating random 16-QAM symbols
Orig_symbols = randi([0, M-1], 1, N_symbols);
Mapped_symbols = qammod(Orig_symbols, M, UnitAveragePower=true);


%% Tapped Delay Line (TDL) Realization

% TASK 2A: GENERATE A TDL REALIZATION (h_TDL)
%
% h_TDL should be a ROW vector with around 3 to 7 values
%
% You can use complex normal distributed channels
% 1/sqrt(2) .* (randn(1, L) + 1i*randn(1, L))
%
% You can also apply some power delay profile with the exp function
%
% You can also just create your own custom vector
%
% Experiment with different different TDL impulse response vectors
h_TDL = 1/sqrt(2) * (randn(1, 3) + 1i*randn(1, 3));


% Normalize so that average signal power is 1
h_TDL = h_TDL/norm(h_TDL);
L = length(h_TDL);


%% Signal Transmission

% Generating noise to match the SNR (average signal power is 1)
sigma_Noise = sqrt(1/SNR);
Noise = sigma_Noise * 1/sqrt(2) .* (randn(1, N_symbols) + 1i*randn(1, N_symbols));

% Transmitting the data through the channel
%
% The filter function filter(b,1,x) applies the finite impulse response (FIR)
% filter defined by b to the vector x (could also be done with conv(b,x))
% 
% Using filter(b,a,x) would result in an infinite impulse response (IIR)
% filter instead
%
Received_symbols = filter(h_TDL, 1, Mapped_symbols) + Noise;


%% Time-Domain Equalization (TDE)

% To estimate the original symbols, we need to somehow reverse the TDL
% echoes. To do this, we simultaneously look at multiple received symbols
% at the same time, and use these to estimate ONE of the original symbols
% that correlates most with the received symbols
% 
% We get a signal model y = Hx + n with the effective channel matrix H

% EXAMPLE
% For a TDL of length 3 ([h_0, h_1, h_2]), and Rx block size 5, H looks like
%
%       h_2  h_1  h_0  0    0    0    0   
%       0    h_2  h_1  h_0  0    0    0
% H =   0    0    h_2  h_1  h_0  0    0
%       0    0    0    h_2  h_1  h_0  0
%       0    0    0    0    h_2  h_1  h_0


% TASK 2B: CREATE THE EFFECTIVE CHANNEL MATRIX
%
% Use the toeplitz function to do this
% 
% H = toeplitz(c,r), where 
% c is the first COLUMN of H and
% r is the first ROW of H (you should use flip(h_TDL) for this).
%
% Remember the zeros! There are (N_r-1) zeros in both c and r 
% (zeros function)

H = toeplitz([h_TDL(end) zeros(1,N_r-1)], [flip(h_TDL) zeros(1,N_r-1)]);


% We can again use either matched filtering (MF), zero forcing (ZF), or
% minimum mean square estimation (MMSE) (similar to MIMO case)


% Matched filter (MF) 
% (calling diag twice returns a matrix with only the diagonal elements left)
% This is basically just H' with some scaling
F_MF = inv(diag(diag(H'*H))) * H';


% Zero forcing (ZF) 
% F is set to be the inverse (or Moore Penrose pseudo-inverse) of H
% (inverting the channel) 
% Note that H is always a wide matrix 
F_ZF = H' * inv(H*H');


% Minimum mean square estimator (MMSE)
N_0 = sigma_Noise^2; % Noise power

F_MMSE = H' * inv(H*H' + N_0*eye(N_r));
% or 
%F_MMSE = inv(H'*H + N_0*eye(N_r+L-1)) * H'; 
 

%% SINR of each Estimated Sample (Lec 6 Appendix)

% Based on a block of received signal samples y_i, we can estimate the
% transmitted signal samples x_i (x_est = F * y)
%
% Each of the samples on x_est will have a different SINR.
% The samples in the middle will have a higher SINR, and the samples on the
% top and bottom have the lowest SINR. 

% Signal and noise powers after filtering
Signal_power_MF = abs(F_MF*H).^2;
Noise_power_MF = N_0*diag(F_MF*F_MF');

Signal_power_ZF = abs(F_ZF*H).^2;
Noise_power_ZF = N_0*diag(F_ZF*F_ZF');

Signal_power_MMSE = abs(F_MMSE*H).^2;
Noise_power_MMSE = N_0*diag(F_MMSE*F_MMSE');

% Loop through the different x_est samples
SINR_vec_MF = zeros(1,N_r+L-1);
SINR_vec_ZF = zeros(1,N_r+L-1);
SINR_vec_MMSE = zeros(1,N_r+L-1);
for k = 1:(N_r+L-1)
    
    % TASK 2C: CALCULATE THE SINR OF THE k:th ESTIMATE
    %
    % Signal_power_FILTER(k,i) gives how much power from the i:th
    % transmitted symbol is received at the k:th estimate for the given filter
    %
    % Noise_power_FILTER(k) gives the noise power for the given filter
    
    % SINR
    SINR_vec_MF(k) = Signal_power_MF(k,k)/(Noise_power_MF(k)+sum(Signal_power_MF(k,:)-Signal_power_MF(k,k)));%POWER_MF/INTERFERENCE_NOISE_MF;
    SINR_vec_ZF(k) = Signal_power_ZF(k,k)/(Noise_power_ZF(k)+sum(Signal_power_ZF(k,:))-Signal_power_ZF(k,k));%POWER_ZF/INTERFERENCE_NOISE_ZF;
    SINR_vec_MMSE(k) = Signal_power_MMSE(k,k)/(Noise_power_MMSE(k)+sum(Signal_power_MMSE(k,:))-Signal_power_MMSE(k,k));%POWER_MMSE/INTERFERENCE_NOISE_MMSE;
end


% SINR plots
figure(5);clf;
hold on
plot(1:(N_r+L-1), 10*log10(SINR_vec_MF))
plot(1:(N_r+L-1), 10*log10(SINR_vec_ZF), '--', 'Linewidth', 2)
plot(1:(N_r+L-1), 10*log10(SINR_vec_MMSE))
legend(["MF", "ZF", "MMSE"], 'Location', 'best')
title("SINR of the estimated signal values")
ylabel("SINR (dB)")
xlabel("Index of the estimated signal value")
xlim([1, (N_r+L-1)])
grid on; grid minor;


%% Choosing the Delay Based on the SINR

% For each block of received signal samples y_i, only one x_j is estimated.
% To choose which one, the highest SINR estimate is used (often in the middle)

% Argmax function (second output of max)
[~, D_MF] = max(SINR_vec_MF);
[~, D_ZF] = max(SINR_vec_ZF);
[~, D_MMSE] = max(SINR_vec_MMSE);


%% Transversal Filter

% Only one element of the estimator vector z = F*y is used
% We can therefore only use the D:th row of F

% TASK 2D: CONSTRUCT A TRANSVERSAL FILTER WHICH CORRESPONDS TO THE HIGHEST
% SINR ESTIMATE

f_vec_MF = F_MF(D_MF, :); % use D_MF
f_vec_ZF = F_ZF(D_ZF, :); % use D_ZF
f_vec_MMSE = F_MMSE(D_MMSE, :); % use D_MMSE


%% Applying the Block Based Time Domain Equalizers

Z_MF = zeros(1,N_symbols-N_r+1);
Z_ZF = zeros(1,N_symbols-N_r+1);
Z_MMSE = zeros(1,N_symbols-N_r+1);

for i=1:N_symbols-N_r+1
    Rx_Block = Received_symbols(i:i+N_r-1).'; 
    Z_MF(i) = f_vec_MF * Rx_Block;
    Z_ZF(i) = f_vec_ZF * Rx_Block;
    Z_MMSE(i) = f_vec_MMSE * Rx_Block;
end

% Ignore the endpoints where the TDL is not fully realized
Z_MF = Z_MF(L:end-L+1);
Z_ZF = Z_ZF(L:end-L+1);
Z_MMSE = Z_MMSE(L:end-L+1);

%% Scatter Plots for the Received and Filtered Symbols

% Received symbols
figure(6); clf;
subplot(2,2,1);
scatter(real(Received_symbols), imag(Received_symbols), 3,'filled')
title('Constellation of Received Symbols')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

% MF filter
subplot(2,2,2);
scatter(real(Z_MF), imag(Z_MF), 3,'filled')
title('Constellation of MF Receiver')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

% ZF filter
subplot(2,2,3);
scatter(real(Z_ZF), imag(Z_ZF), 3,'filled')
title('Constellation of ZF Receiver')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

% MMSE filter
subplot(2,2,4);
scatter(real(Z_MMSE), imag(Z_MMSE), 3,'filled')
title('Constellation of MMSE Receiver')
xlabel("real"); ylabel("imag");
grid on; grid minor; axis equal; axis square;

%% Symbol Error Rates

% MF
Symbol_estimate_MF = qamdemod(Z_MF, M, UnitAveragePower=true);
Symbol_error_rate_MF = 1 - mean(Symbol_estimate_MF == Orig_symbols(D_MF:D_MF+length(Symbol_estimate_MF)-1))

% ZF
Symbol_estimate_ZF = qamdemod(Z_ZF, M, UnitAveragePower=true);
Symbol_error_rate_ZF = 1 - mean(Symbol_estimate_ZF == Orig_symbols(D_ZF:D_ZF+length(Symbol_estimate_ZF)-1))

% MMSE
Symbol_estimate_MMSE = qamdemod(Z_MMSE, M, UnitAveragePower=true);
Symbol_error_rate_MMSE = 1 - mean(Symbol_estimate_MMSE == Orig_symbols(D_MMSE:D_MMSE+length(Symbol_estimate_MMSE)-1))














