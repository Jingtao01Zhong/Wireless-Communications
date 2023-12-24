%% Task 2 (Open-Loop and Closed-Loop MIMO Capacity)


%%
clear all;


%% Initialization (Feel free to experiment with these values)

% Number of Tx (transmitter) antennas
N_t = 8;

% Number of Rx (receiver) antennas
N_r = 8;


% SNRs to evaluate
SNR_dB_vec = -20:5:20;
SNR_vec=10.^(SNR_dB_vec./10);


% Number of channel realizations for Monte Carlo
N_channel_realizations = 200;


%% CHANNEL REALIZATIONS 

% H_array(:,:,j) gives the j:th N_r x N_t channel matrix realization
H_array = 1/sqrt(2) .* (randn(N_r, N_t, N_channel_realizations) + 1i*randn(N_r, N_t, N_channel_realizations));

% Average energy of received signal
% This is given by the expected value of the diagonal values of H*H', which is 
% just the number of transmitter antennas (as P_t = 1)
E_s = N_t;

% Average transmission power is 1
P_t = 1;

%% Ergodic Capacity of Rayleigh Fading MIMO Channnels Using Monte Carlo

% Keep track of the average/ergodic capacities and the variance for each SNR
% OPEN LOOP
Open_loop_Capa_avg = zeros(size(SNR_vec));


% CLOSED LOOP (WITH power allocation)
Closed_loop_Capa_avg = zeros(size(SNR_vec));


% Loop through the different SNRs
for k = 1:length(SNR_vec)
    
    SNR = SNR_vec(k);
 
    % Noise power
    N_0 = E_s/SNR;
    

    % Keep track of the capacities for each realization (MC = Monte Carlo)
    OL_Capa_MC = zeros(N_channel_realizations, 1);
    CL_Capa_MC = zeros(N_channel_realizations, 1);

    % Loop through the different channel realizations
    for j = 1:N_channel_realizations
        
        % Take the j:th channel realization
        H = H_array(:,:,j);
        % The (i,j) element of H is the channel between the 
        % i:th RECEIVER antenna and the j:th TRANSMITTER antenna (y = Hx + n).
        %
        % Each of these channels is assumed to be an independent Rayleigh channel.


        % SINGULAR VALUE DECOMPOSITION (SVD)
        % The svd function returns two orthogonal matrices U,V 
        % and a diagonal matrix Sigma such that U*Sigma*V' = H
        [U, Sigma, V] = svd(H);
        
        
        % OPEN-LOOP MIMO
        
        % Perfect channel state information (CSI) at the receiver (but not the transmitter)
        %
        % The receiver can filter the signal with the unitary matrix U'. 
        % Note that noise characteristics are not affected, as U is unitary
        %
        % The effective channel is thus
        % H_open_loop = U'*H = U'*(U*Sigma*V') = (U'*U)*Sigma*V' = Sigma*V'
        %
        % There are r = min(N_t, N_r) channels, with the gain of 
        % each channel given by the singular values of H (diagonal of Sigma).
        %
        % The transmitter has no idea about the channel matrix, so the best
        % it can do is to divide its total power evenly accross all
        % antennas. The transmitter is forced to use power on all
        % "directions" V, even the ones that don't reach the transmitter
        % (or the ones that correspond to small singular values).

        H_open_loop = Sigma*V'; %(= U'*H)

        % TASK 2.1 Calculate the open-loop MIMO channel capacity using the 
        % formula from Lecture 9, slide 20. 
        % You can either calculate the capacity using the matrix H (or H_open_loop)
        % directly, or use the singular values from the diagonal of Sigma.
        %
        % (total transmission power P_t = 1, noise power is given by N_0,
        % and number of transmitter antennas is given by N_t)
        C_subchannel = log2(diag(ones(N_t,1))+P_t/N_t/N_0*Sigma.*Sigma);
        OL_Capa_MC(j) = trace(C_subchannel);
        
        

        % CLOSED-LOOP MIMO
        
        % Perfect channel state information (CSI) at the transmitter and the receiver.
        %
        % The transmitter can thus precode the signal with the unitary matrix V from the svd.
        % The receiver can again filter the signal with the unitary matrix U'. 
        % Note that noise characteristics are not affected, as U and V are both unitary 
        %
        % The effective channel is now just a diagonal matrix of the singular values of H
        % H_closed_loop = U'*H*V = U'*(U*Sigma*V')*V = (U'*U)*Sigma*(V'*V) = Sigma
        %
        % There are r = min(N_t, N_r) channels, with the gain of 
        % each channel given by the singular values of H (diagonal of Sigma).
        
        H_closed_loop = Sigma; %(= U'*H*V)
        
        % This means that the transmitted symbols have no 
        % Intersymbol Interference (ISI) at the receiver.
        % Thus an optimal receiver is very easy to implement.

        % However, if power is allocated evenly, the theoretical capacity
        % is EQUAL to the open loop case. To gain capacity, we need to 
        % allocate power between the channels.
        
        
        % WATERFILLING ALGORITHM (see Lecture 10, slide 14)
        
        % To allocate the power, we use the waterfilling algorithm.

        % Normalized gain
        g = diag(Sigma).^2/N_0;
        % For the algorithm, we need the inverse values of the gains
        inv_g = 1./g;

        % Number of parallel channels
        r = length(g);
        % Available power (=1)
        P_t = 1;
        
        % For the power, we have the constraints 
        % P_i >= 0 and sum_i (P_i) = P_t
        
        % The optimal power allocation is obtained using Lagrangian optimization
        % and the so called KKT conditions.
        
        % Keep track of the channels that are used (all at the start)
        used_channels = true(r,1);
        P_alloc = zeros(r, 1);
        
        % Iterate until a suitable power allocation is found
        while 1
            % ONLY ALLOCATE POWER TO USED CHANNELS
            N_used_channels = sum(used_channels);
    
            % TASK 2.2 Calculate 1/lambda (1/λ) for the waterfilling algorithm.
            % Use inv_g(used_channels) to get the inverse gains of the
            % used channels.
            % Remember to use N_used_channels instead of N_t.
            one_per_lambda = (P_t + sum(inv_g(used_channels)))/N_used_channels;
            
            % TASK 2.3 Calculate the power allocation for the used channels 
            % using 1/λ. Remember to only update the used channels.
            P_alloc(used_channels) = one_per_lambda - inv_g(used_channels);
            

            % If some channel gets negative power, drop the channel and set its
            % power to zero, and reiterate without the channel
            negative_channels = P_alloc < 0;
            if sum(negative_channels)==0
                % No negative channel powers, valid power allocation
                break % (gets out of the while loop)
            else
                % We need to drop the negative channels and reiterate
                P_alloc(negative_channels) = 0;
                used_channels(negative_channels) = 0;
            end
        end
        
        % We now have a power allocation (P_alloc) that maximizes the
        % total capacity across channels.
        
        % TASK 2.4 Calculate the SNR of each parallel channel using the 
        % power allocation (P_alloc) and the normalized gain (g).
        SNR_P_alloc = P_alloc .* g;
        
        
        % TASK 2.5 Use the SNRs of each parallel channel to calculate the 
        % capacity of the closed loop MIMO channel by using Shannon's 
        % formula (C = log2(1+SNR)) separately for each parallel channel.
        CL_Capa_MC(j) = sum(log2(1+SNR_P_alloc));

    end

    % We can now calculate the average/ergodic capacity
    % OPEN LOOP
    Open_loop_Capa_avg(k) = mean(OL_Capa_MC);
    
    
    % CLOSED LOOP (WITH power allocation)
    Closed_loop_Capa_avg(k) = mean(CL_Capa_MC);
end



%% Plotting the Capacities

% Plotting capacity
figure(4); clf
subplot(2,1,1)
hold on
plot(SNR_dB_vec, Open_loop_Capa_avg)
plot(SNR_dB_vec, Closed_loop_Capa_avg)
grid on
grid minor
xlabel('SNR (dB)')
ylabel('Ergodic Capacity (bits per channel use)')
title('Open-Loop vs. Closed-Loop MIMO N_t = '+string(N_t)+', N_r = '+string(N_r))
legend(["Open-Loop MIMO", "Closed-Loop MIMO (with Power Allocation)"],'Location', 'best')
xlim([min(SNR_dB_vec),max(SNR_dB_vec)])

% Plotting the ratio between capacities 
ratio = Closed_loop_Capa_avg./Open_loop_Capa_avg;
subplot(2,1,2)
hold on
plot(SNR_dB_vec, ratio)
grid on
grid minor
xlabel('SNR (dB)')
ylabel('Ergodic Capacity Ratio')
legend("$\frac{C_{Closed}}{C_{Open}}$",'Location', 'best', 'Interpreter','Latex','fontsize',20)
xlim([min(SNR_dB_vec),max(SNR_dB_vec)])


%% Random Matrix Singular Values

% Number of singular values/channels
r = min(N_r, N_t);

% Keep track of the singular values
Sigma_average = zeros(r,1);
% Keep track of the average power of each singular value channel
Sigma_squared_average = zeros(r,1);

% Loop through the different channel realizations
for j = 1:N_channel_realizations
    
    % Take the j:th channel realization
    H = H_array(:,:,j);

    % SINGULAR VALUE DECOMPOSITION (SVD)
    % This time we only want the singular values
    % These are ordered from largest to smallest
    sigma = svd(H);

    % Keep track of the sums
    Sigma_average = Sigma_average + sigma;
    Sigma_squared_average = Sigma_squared_average + sigma.^2;
end

% Divide by N_channel_realizations to get the averages
Sigma_average = Sigma_average/N_channel_realizations;
Sigma_squared_average = Sigma_squared_average/N_channel_realizations;

% Pad the vector with N_t-r zeros to include the transmit directions that
% don't reach the receiver
Sigma_average = [Sigma_average; zeros(N_t-r,1)];
Sigma_squared_average = [Sigma_squared_average; zeros(N_t-r,1)];

figure(5); clf
hold on
plot(1:N_t, Sigma_average)
plot(1:N_t, Sigma_squared_average)
grid on
grid minor
xlabel('Singular value index')
title('Random Matrix Average Singular Values')
legend(["Average of n:th singular value", "Average of n:th squared singular value"],'Location', 'best')
xlim([1, N_t])

















