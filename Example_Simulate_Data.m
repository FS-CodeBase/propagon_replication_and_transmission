% SCRIPT: EXAMPLE_SCRIPT
% AUTHOR: Fabian Santiago
% E-mail: fabiansantiago707@gmail.com
% DATE: 11/15/2021
% DESCRIPTION: This script shows how the functions are used to generate 
%              simulated data and plots the generated data.


% Add code folder to Path
addpath('./code/') % Model solutions, likelihood, recursive log of sum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                Generate Simulated Data               %%%%%%%%%%

lambda = 0.70; % Prion replication rate
rho    = 0.30; % Transmission bias
sampling_times = 0:7; % Sample every hour for 7 hours after start
sampling_rate  = 32;  % 32 samples per hour

% Generate Simulated Data and save in simulated_data folder.
% Data is also saved as:
% simulated_data/simdata_mu10p0sig1p0lambda0p70rho0p30smph32p0.mat
[propagon_data,sampling_times] = ...
        simulate_data(lambda,rho,sampling_times,sampling_rate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                  Plot Simulated Data                 %%%%%%%%%%
% use latex interpreter
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

plot_one_data(propagon_data,sampling_times)
title(['$\lambda = ',num2str(lambda),'$ and $\rho = ',num2str(rho),'$'],...
        'Interpreter','latex','FontSize',20)