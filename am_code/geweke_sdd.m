function p_val = geweke_sdd(in_signal)
% FUNCTION: GEWEKE_SDD computes Geweke's spectral density diagnostic.
% % AUTHOR: Fabian Santiago
% % EMAIL: fsantiago3@ucmerced.edu
% % DATE: 11/23/2020
% % INPUTS: IN_SIGNAL an r by c matrix where each c is one parameter and r
% %                   is an mcmc sample.

% Determine N, N_a, N_b, and N_star (determined randomly)

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