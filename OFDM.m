%% Task 1 (OFDM)

%% Initialization (Feel free to experiment with these values)

% Length of tapped delay line (TDL)
L = 3;
% Decay parameter for the TDL (should be positive) (lambda = 0 means no decay)
% (Higher value = faster decay)
lambda = 0.5;

% OFDM parameters
N_Carrier = 64; % Number of carriers (size of fft) (Should be a power of 2)

N_CP = L-1; % Length of syclic prefix (needs to be L-1)
Frame_efficiency = N_Carrier/(N_Carrier+N_CP);

% SNRs to evaluate
SNR_dB_vec = 0:10:70;
SNR_vec = 10.^(SNR_dB_vec./10);

% Number of OFDM symbols/blocks
N_OFDM_symbols = 1e4;

% Number of symbols
N_symbols = N_OFDM_symbols * N_Carrier;

%% Gray constellation for M-QAM 

% Switch between different QAM symbols (M should be a power of 2)
%M = 2 % BPSK
%M = 4 % QPSK
M = 16; % 16-QAM

% Constellation scatter plot
Constellation_QAM = qammod(0:M-1, M, UnitAveragePower=true);

figure(1); clf;
hold on
scatter(real(Constellation_QAM(:)), imag(Constellation_QAM(:)), 30,'filled')
line([-1.5 1.5], [0 0], 'Color', 'black')
line([0 0], [-1.5 1.5], 'Color', 'black')
xlabel("real"); ylabel("imag");
title(string(M)+'-QAM constellation')
grid on; grid minor; axis square;


%% Data Deneration

% Generating random M-QAM symbols
Original_symbols = randi([0, M-1], 1, N_symbols);
Mapped_symbols = qammod(Original_symbols, M, UnitAveragePower=true);


%% Data Transmission

% Keep track of the equalized symbols for each SNR
Equalized_symbols = zeros(length(SNR_vec), N_symbols);
for j = 1:length(SNR_vec)
    SNR = SNR_vec(j);
    sigma_Noise = sqrt(1/SNR);
    % Loop through the different OFDM symbols/blocks
    for k = 1:N_OFDM_symbols
        
        % Tapped Delay Line (TDL) Realization 
        % (Rayleigh fading channel taps with decreasing intensity)
        h_TDL = 1/sqrt(2) .* (randn(1, L) + 1i*randn(1, L)) .* sqrt(exp(-lambda.*(1:L)));
        % Normalize so that the average signal power is 1
        h_TDL = h_TDL/norm(h_TDL);

        % TRANSMITTER PROCESSING

        % 1) Data block to be sent (frequency domain data)
        x_freq = Mapped_symbols(1+(k-1)*N_Carrier : k*N_Carrier);

        % 2) Transferring the frequency domain data into time domain data

        % TASK 1.1: Use the ifft-function to transfer the frequency domain 
        % data into time domain data. Remember to scale by sqrt(N_Carrier)
        % so that the average transmit power stays at 1.
        x_time = sqrt(N_Carrier)*ifft(x_freq);

        % 3) Adding the cyclic prefix
        % TASK 1.2: Add a cyclic prefix of length (L-1) to make the
        % effective channel matrix circular
        prefix = x_time(N_Carrier-L+2:end);
        OFDM_symbol = [prefix, x_time];

       
        % 4) Transmission through the Tapped Delay Line
        Noise = sigma_Noise * 1/sqrt(2) .* (randn(1, N_CP+N_Carrier) + 1i*randn(1, N_CP+N_Carrier));
        Received_OFDM_symbol = filter(h_TDL,1,OFDM_symbol) + Noise;

        
        % RECEIVER PROCESSING

        % 1) Cyclic prefix removal
        % TASK 1.3: Remove the cyclic prefix from the received OFDM symbol.
        Received_OFDM_symbol = Received_OFDM_symbol(L:end);

        % 2) Transferring the data back into the frequency domain
        % TASK 1.4: Use the fft-function to transfer the time domain 
        % data back into the frequency domain. 
        % This time the result needs to be scaled by sqrt(1/N_Carrier).
        y_freq = 1/sqrt(N_Carrier).*fft(Received_OFDM_symbol);
    
        % 3) Calculate the frequency domain channel matrix elements
        % TASK 1.5: Use the fft-function to calculate the elements of the
        % diagonalized channel matrix.
        % You need to use the variables h_TDL and N_Carrier for this.
        h_diag = fft(h_TDL,N_Carrier);

        % 4) Apply equalization in the frequency domain
        % TASK 1.6: Use h_diag to equalize the channel in the frequency
        % domain. You should use ZF for this. 
        % (We essentially have y_freq = h_diag .* x_freq + noise, 
        % where x_freq is the thing we want to predict.)
        Equalized_block = h_diag.^-1 .* y_freq;

        % 5) Store the equalized symbols
        Equalized_symbols(j, 1+(k-1)*N_Carrier : k*N_Carrier) = Equalized_block;
    end
end

%% Symbol Error Rates

% Demodulation
Symbol_estimate = qamdemod(Equalized_symbols, M, UnitAveragePower=true);

% Symbol error for each SNR
Symbol_error_rate = 1 - mean(Original_symbols == Symbol_estimate,2);

% Plotting symbol error against SNR
figure(2); clf;
hold on
plot(SNR_dB_vec, Symbol_error_rate)
set(gca, 'YScale', 'log')
grid on
grid minor
xlabel('SNR (dB)')
ylabel('SER')
title('Symbol Error Rate (SER) using OFDM with '+string(M)+'-QAM symbols.')
xlim([min(SNR_dB_vec), max(SNR_dB_vec)])






















