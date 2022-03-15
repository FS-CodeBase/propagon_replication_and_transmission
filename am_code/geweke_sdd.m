function p_val = geweke_sdd(in_signal)
% FUNCTION: GEWEKE_SDD computes Geweke's spectral density diagnostic.
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
% 
% DESCRIPTION
% % Computes Gewekes convergence diagnostic 
%
% INPUTS 
% % in_signal: an r by c matrix where each c is one parameter and r
% %            is an mcmc sample.
%
% OUTPUTS
% % p-value associated with Geweke's convergence diagnostic value

% Determine N, N_a, N_b, and N_star
N = size(in_signal,1);
N_A = ceil(N/10); N_B = ceil(N/2); N_str = N-N_B+1; 

% Estimate the spectral density
S_A = periodogram(in_signal(1:N_A));
S_B = periodogram(in_signal((N_str:N)));

% Estimate the mean 
theta_A = mean(in_signal(1:N_A)); 
theta_B = (N_str/N_B)*mean(in_signal(N_str:N));

% Z-statistic estimate
Zstat = (theta_A-theta_B)/sqrt(1/N_A*S_A(1)+1/N_B*S_B(1));
p_val = 1-2*normcdf(-abs(Zstat),0,1);