%% Task 1 (AWGN Channel)

% In this task, you evaluate the effect of Forward Error Correction (FEC) 
% on the Bit Error Rate (BER) of an AWGn channel.

%%
clear all;

%% Initialization

% Number of bits to be sent (should be at least 1e3)
N_bits = 1e5;

% Initializing the QPSK constellation
Modulation = 'QPSK';
QPSK_constellation = [1+1i, -1+1i, 1-1i, -1-1i]/sqrt(2);
Bits_per_symbol = log2(length(QPSK_constellation)); %(=2)

% SNRs to evaluate (feel free to change)
SNRdB_vec = -2:10;


%% Bit generation

% Generating random data
Data_bits = randi([0 1], 1,N_bits);


%% 1) No Coding, Data Sent Directly

% The data bits are immediately mapped to QPSK symbols
% Reshaping data bits so that each QPSK symbol correspons to 2 bits
N_QPSK = length(Data_bits)/2;
Reshaped_bits = reshape(Data_bits, 2, N_QPSK);

% Mapping the bits to QPSK symbols
QPSK_symbols = QPSK_constellation([1 2]*Reshaped_bits+1);

% We keep track of the bit error rate (BER)
BER_no_coding = zeros(1, length(SNRdB_vec));

% Transmitting the QPSK symbols through AWGN channel with varying SNR
for k = 1:length(SNRdB_vec)
    
    % Computing Power of Thermal Noise (Average signal power is set to 1)
    SNRdB = SNRdB_vec(k);
    SNR = 10^(SNRdB/10);
    
    % STEP 1
    % Generate complex gaussian noise with the specified noise power
    Noise_Power=1/SNR; 
    Noise = sqrt(Noise_Power) * sqrt(0.5)*(randn(1,N_QPSK)+1j*randn(1,N_QPSK)); % FILL THIS

    % Transmitting symbols (adding noise)
    Received_QPSK_symbols = QPSK_symbols + Noise;

    % Demodulation
    Received_bits = 1*[real(Received_QPSK_symbols)<0;
                       imag(Received_QPSK_symbols)<0];
    Received_bits = reshape(Received_bits, 1, []);

    % There is no decoding 
    Decoded_bits = Received_bits;
    
    % STEP 2
    % Calculate the bit error rate (BER) by comparing the
    % vectors Decoded_bits and Data_bits
    BER_no_coding(k) = sum(Decoded_bits ~= Data_bits)./length(Decoded_bits); % FILL THIS

end

BER_no_coding



%% 2) Repetition Coding

% Bits are encoded by repeating them r (=3) times
% Repetition amount
r = 3;

% Coding rate
CodeRate_Repetition = 1/r; 

% Encoding (repeating each bit r times)
Encoded_bits = kron(Data_bits, ones(1,r));

% Reshaping the encoded bits so that each QPSK symbol correspons to 2 bits
N_QPSK = length(Encoded_bits)/2;
Reshaped_bits = reshape(Encoded_bits, 2, N_QPSK);

% Mapping the bits to QPSK symbols
QPSK_symbols = QPSK_constellation([1,2]*Reshaped_bits+1);

% We keep track of the bit error rate (BER) and throughput
BER_repetition_coding = zeros(1, length(SNRdB_vec));
Throughput_repetition_coding = zeros(1, length(SNRdB_vec));

% Transmitting the QPSK symbols through AWGN channel with varying SNR
for k = 1:length(SNRdB_vec)
    
    % Computing Power of Thermal Noise (Average signal power is set to 1)
    SNRdB = SNRdB_vec(k);
    SNR = 10^(SNRdB/10);

    % STEP 1
    % Generate complex gaussian noise with the specified noise power
    Noise_Power=1/SNR; 
    Noise = sqrt(Noise_Power) * sqrt(0.5)*(randn(1,N_QPSK)+1j*randn(1,N_QPSK)); % FILL THIS

    % Transmitting symbols (adding noise)
    Received_QPSK_symbols = QPSK_symbols + Noise;

    % Demodulation
    Received_soft_bits = [real(Received_QPSK_symbols)
                          imag(Received_QPSK_symbols)];
    Received_soft_bits = reshape(Received_soft_bits, 1, []);

    % These soft bits now have to be decoded
    % Group the estimated bits to groups of r
    Grouped_soft_bits = reshape(Received_soft_bits, r, []);

    % Summing the grouped bits and seeing if the total is 
    % positive (bit = 0) or negative (bit = 1)
    Decoded_bits = sum(Grouped_soft_bits, 1) < 0;

    % STEP 2
    % Calculate the bit error rate (BER) by comparing the
    % vectors Decoded_bits and Data_bits
    BER_repetition_coding(k) = sum(Decoded_bits ~= Data_bits)./length(Decoded_bits); 
    % FILL THIS

    % STEP 3
    % Use this BER to compute the throughput for this SNR
    % (Use the values CodeRate_Repetition (=1/3) and Bits_per_symbol (=2))
    Throughput_repetition_coding(k) = CodeRate_Repetition .* Bits_per_symbol .* (1-BER_repetition_coding(k)); %FILL THIS
end

BER_repetition_coding 



%% 3) Turbo Coding

% Setting up the turbo encoder/decoder pair
frameLen = 1000; % Frame length
N_frames = N_bits/frameLen;
intrlvrIndices = randperm(frameLen);
NumIterations = 4;
Turbo_Encoder = comm.TurboEncoder('InterleaverIndices',intrlvrIndices);
Turbo_Decoder = comm.TurboDecoder('InterleaverIndices',intrlvrIndices,'NumIterations',NumIterations);

% The codeRate is 1/3 for this code
CodeRate_Turbo = 1/3;

% Encoding the bits 
Encoded_bits = zeros(1, 3*N_bits);
for i=1:N_frames
    Data_frame = Data_bits(1+(i-1)*frameLen:i*frameLen)'; 
    Encoded_frame = step(Turbo_Encoder, Data_frame)';
    Encoded_bits(1+3*(i-1)*frameLen:3*i*frameLen) = Encoded_frame(1:end-12);
end

%% Transmitting the Encoded Data

% Reshaping the encoded bits so that each QPSK symbol correspons to 2 bits
N_QPSK = length(Encoded_bits)/2;
Reshaped_bits = reshape(Encoded_bits, 2, N_QPSK);

% Mapping the bits to QPSK symbols
QPSK_symbols = QPSK_constellation([1,2]*Reshaped_bits+1);

% We keep track of the bit error rate (BER) and throughput
BER_turbo_coding = zeros(1, length(SNRdB_vec));
Throughput_turbo_coding = zeros(1, length(SNRdB_vec));

% Transmitting the QPSK symbols through AWGN channel with varying SNR
for k = 1:length(SNRdB_vec)
    % Computing Power of Thermal Noise (Average signal power is set to 1)
    SNRdB = SNRdB_vec(k);
    SNR = 10^(SNRdB/10);

    % STEP 1
    % Generate complex gaussian noise with the specified noise power
    Noise_Power=1/SNR; 
    Noise = sqrt(Noise_Power) * sqrt(0.5)*(randn(1,N_QPSK)+1j*randn(1,N_QPSK)); % FILL THIS

    % Transmitting symbols (adding noise)
    Received_QPSK_symbols = QPSK_symbols + Noise;
    

    % Demodulation
    Received_soft_bits = [real(Received_QPSK_symbols)
                          imag(Received_QPSK_symbols)];
    Received_soft_bits = reshape(Received_soft_bits, 1, []);

    % These soft bits now have to be decoded
    % This is done using the Turbo decoder
    
    % The Turbo decoder uses log-likelihood ratios (LLR)
    % ln(P(bit=1)/P(bit=0)) to do the decoding
    noisePower_per_branch = Noise_Power/2; % Evenly divided between the real and imaginary parts
    LLR = (1/(2*noisePower_per_branch)) .* ((Received_soft_bits - sqrt(1/2)).^2 - (Received_soft_bits + sqrt(1/2)).^2);

    Decoded_bits = zeros(1, N_bits);
    for i=1:N_frames
        LLR_frame = [LLR(1+3*(i-1)*frameLen:3*i*frameLen)'; zeros(12,1)]; 
        Decoded_frame = step(Turbo_Decoder,LLR_frame)';
        Decoded_bits(1+(i-1)*frameLen:i*frameLen) = Decoded_frame;
    end

    % STEP 2
    % Calculate the bit error rate (BER) by comparing the
    % vectors Decoded_bits and Data_bits
    BER_turbo_coding(k) = sum(Decoded_bits ~= Data_bits)./length(Decoded_bits);

    % STEP 3
    % Use this BER to compute the throughput for this SNR
    % (Use the values CodeRate_Turbo (=1/3) and Bits_per_symbol (=2))
    Throughput_turbo_coding(k) = CodeRate_Repetition .* Bits_per_symbol .* (1-BER_turbo_coding(k));
end

BER_turbo_coding

 


%% Plotting the Bit Error Rates (BER) 
figure(1); clf;
hold on
plot(SNRdB_vec, BER_no_coding)
plot(SNRdB_vec, BER_repetition_coding)
plot(SNRdB_vec, BER_turbo_coding)
set(gca, 'YScale', 'log')
grid on
grid minor
xlabel('SNR (dB)')
ylabel('BER')
title('Bit error rate (BER) using QPSK in an AWGN channel')
legend(["No encoding", "Repetition code (r = "+ string(r)+ ")", "Turbo code (Code rate 1/3)"])
xlim([min(SNRdB_vec), max(SNRdB_vec)])


%% Capacity of the Channel

Channel_capacity_vec = zeros(1, length(SNRdB_vec));

for k = 1:length(SNRdB_vec)
    % Converting from dB to linear
    SNRdB = SNRdB_vec(k);
    SNR = 10^(SNRdB/10);

    % STEP 4
    % Compute the capacity for the AWGN channel using Shannon's formula (Lec 1, Slide 6)
    Channel_capacity_vec(k) = log2(1 + SNR)% FILL THIS
end
Channel_capacity_vec


%% Plotting Throughput and Capacity
figure(2); clf
hold on
plot(SNRdB_vec, Throughput_repetition_coding)
plot(SNRdB_vec, Throughput_turbo_coding)
plot(SNRdB_vec, Channel_capacity_vec)
grid on
grid minor
xlabel('SNR (dB)')
ylabel('Throughput (bits per transmission)')
title('Throughput for repetition and turbo code in an AWGN channel')
legend(["Repetition code (r = "+ string(r)+ ")", "Turbo code (Code rate 1/3)", "AWGN channel capacity"],'Location', 'SE')
xlim([min(SNRdB_vec),max(SNRdB_vec)])
ylim([0,1.2])




