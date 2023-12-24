% In this Matlab example we want to check the Central limit theorem (CLT). 
% Consider the transmitted signal to be a single tone with frequency 1, i.e., sin(2*pi*t). 
% At the receiver, we will receive the transmitted signal from different directions with
% different delay. R = sum(A(xi)*(cos(2*pi*phi_i)+1j*sin((2*pi*phi_i)))).
function ChannelRealization_Task1
%% Settings

clear; close all; clc;
Ndrops = 1e4; % Number of different channel realizations
SampN  = 1e4; % Number of different directions/amplitudes
% For OK experimental distributions, you need at least 1000, for good, 10000 


%% Channel Realization Generation

% Uniformly distributed amplitudes and phases
Amp   = rand(SampN, Ndrops);
Phase = rand(SampN, Ndrops);

% Calculating the contributions of different paths
% (exp(i*x) = cos(x) + i * sin(x))
R = Amp.*(cos(2*pi*Phase) + 1j*sin(2*pi*Phase));


% Summing the paths and setting the expected mean power (|R|^2) to 1
R = sum(R,1)*sqrt(3/SampN);
mean_channel_power = mean(abs(R).^2)


%% Comparing the Real and Imaginary Parts of the Signal Samples to the Normal Distribution

% Observe that the average power is equally divided between the
% real and imaginary parts (0.5 for each)

% Plotting the real parts of the signal samples
figure(1); clf;
subplot(2,1,1);
Nbins=round(Ndrops/10);
histogram(real(R),Nbins,'Normalization','pdf')
hold on, grid on ,grid minor

% Plotting the PDF of the Normal distribution
NormalPDF= @(x) 1/sqrt(2*pi*0.5)*exp(-x.^2);
fplt=fplot(NormalPDF);
fplt.LineWidth=2;
ylim([0 1])
xlim([-3 3])
ylabel('Probability Density')
xlabel('real(R)')
title(['Real Parts of Signal Samples' newline 'Compared to the Normal Distribution'])
legend('Real Parts of Signal Samples','Normal Distribution; \mu=0, \sigma^2=0.5')


% Plotting the imaginary parts of the signal samples
subplot(2,1,2)
histogram(imag(R),Nbins,'Normalization','pdf')
hold on, grid on ,grid minor

% Plotting the PDF of the Normal distribution
fplt = fplot(NormalPDF);
fplt.LineWidth = 2;
ylim([0 1])
xlim([-3 3])
ylabel('Probability Density')
xlabel('im(R)')
title(['Imaginary Parts of Signal Samples' newline 'Compared to the Normal Distribution'])
legend('Imaginary Parts of Signal Samples','Normal Distribution; \mu=0, \sigma^2=0.5')


%% Scatter Plot of Signal Samples

figure(2); clf;
scatter(real(R),imag(R),3,'filled')
axis equal
xlim([-3 3])
ylim([-3 3])
xlabel('real(R)')
ylabel('im(R)')
title('Scatter of Signal Samples')
grid on, grid minor


%% Comparing the Signal Sample Powers to the Exponential Distribution

% Plotting the signal sample powers
figure(3); clf;
subplot(2,1,1)
histogram(abs(R).^2,Nbins,'Normalization','pdf')
hold on, grid on ,grid minor

% Defining the PDF of the exponential distribution (The mean power is 1)
ExponentialPDF = @(x) exp(-x); %FILL THIS


% Plotting the PDF of the exponential distribution
fplt=fplot(ExponentialPDF);
fplt.LineWidth=2;
ylim([0 1])
xlim([0 3])
title(['Signal Sample Powers' newline 'Compared to the Exponential Distribution'])
legend('Signal Sample Powers','Exponential Distribution; \lambda=1')


%% Comparing the Signal Sample Amplitudes to the Rayleigh Distribution.

% Plotting the signal sample amplitudes
subplot(2,1,2)
histogram(abs(R),Nbins,'Normalization','pdf')
hold on, grid on ,grid minor

% Defining the PDF of the Rayleigh distribution (The mean power is 1)
RayleighPDF = @(x) 2*x*exp(-(x^2)); %FILL THIS


fplt = fplot(RayleighPDF);
fplt.LineWidth=2;
ylim([0 1])
xlim([0 3])
title(['Signal Sample Amplitudes' newline 'Compared to the Rayleigh Distribution'])
legend('Signal Sample Powers','Rayleigh Distribution; \sigma^2=1/2')

end