% In this Matlab example, we want to study delay spread, coherence
% bandwidth, and frequency domain autocorrelation.
function DelaySpreadBCCalculation_Task2
%% Settings

clear; close all; clc;
DelayProfParam = 7; %should be greater than 2. Higher value = slower decay


%% Power Delay Profile (PDP)

% Generating and plotting a PDP that decays exponentially
Ntap  = 256; 
Delay = 0:Ntap-1;
PDP   = exp(-Delay/DelayProfParam);
PDP   = PDP/sum(PDP); % normalization
figure(1); clf;
plot(Delay, PDP);
xlim([0,50]);
title('Power Delay Profile Example');
xlabel('Delay')
ylabel('Amplitude Gain')
grid on, grid minor


%% Calculating the RMS delay Spread and the Coherence Bandwith

% Average delay spread
AvgDS = sum(Delay.*PDP);

%RMS Delay spread (FILL THIS)
rmsDS = sqrt(sum(((Delay-AvgDS).^2) .* PDP));


%Coherence bandwidth (FILL THIS)
BC = 1/rmsDS;


% Plotting the autocorrelation in frequency (Fourier transform of PDP)
figure(2); clf; 
FreqResp = fftshift(fft(PDP));
freq = (Delay/Ntap-0.5);
plot(freq, abs(FreqResp))
xlim([-0.5,0.5])
grid on, grid minor
hold on

% Plotting the Coherence Bandwith
[~, BC_loc] = min(abs(freq-BC));
line([BC, BC],[0,abs(FreqResp(BC_loc))], 'Color','red')
line([-BC, -BC],[0,abs(FreqResp(BC_loc))], 'Color','red')
line([-BC, BC],[abs(FreqResp(BC_loc)), abs(FreqResp(BC_loc))], 'Color','red')
legend(["Autocorrelation in Frequence Space", "Coherence bandwith"], 'location', 'SouthEast')
title('Coherence Bandwith in the Generated PDP')
xlabel('Frequency [1/T] (T = delay between taps)')
ylabel('Autocorrelation')


%% Generating a Channel Based on the PDP

% Generating normally distributed complex path gains (mean power = 1)
h = sqrt(0.5)*(randn(1,Ntap)+1j*randn(1,Ntap));

% Change h so that the power of each tap corresponds to the generated PDP
h = h .* sqrt(PDP); % (FILL THIS)


% Plotting the channel in frequency space
figure(3); clf;
Fh = fftshift(fft(h));
freq = (Delay/Ntap-0.5);
plot(freq, abs(Fh).^2)
xlim([-0.5,0.5])
xlabel('Frequency [1/T] (T = delay between taps)')
ylabel('Channel Frequency Response Power')
title('Channel Frequency Response')
grid on, grid minor


end