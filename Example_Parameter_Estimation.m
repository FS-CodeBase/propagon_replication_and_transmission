% SCRIPT: EXAMPLE_PARAMETER_ESTIMATION
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
%
% DESCRIPTION
% % This script shows how the functions are used to perform parameter 
% % estimation using the simulated data and secondly using the 
% % experimental data.
%

% Add adaptive metropolis code to path
addpath('./code/'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         Parameter Estimation: Simulated Data         %%%%%%%%%%
% Load dataset: Weak
load('simulated_data/simdata_mu10p0sig1p0lambda0p70rho0p30smph32p0.mat',...
                                         'propagon_data','sampling_times')

% Setup initial parameter values to begin parameter estimation.
Prms0 = [0.1,0.5]; % Initial parameter estimates [lambda, rho]
Cn    = sqrt([0.1406,0;0,0.0156]); % Initial covariance matrix
num_itrs = 6*10^3;

% Perform parameter estimation using Adaptive Metropolis Algorithm
[PRMS,~] = am_alg(propagon_data,sampling_times,Prms0,Cn,num_itrs);

% Plot adaptive Metropolis chain
figure
T = tiledlayout(2,1);
nexttile
plot(PRMS(:,1),'LineWidth',2,'Color',[0,0,.75])
grid on
legend(['$\bar\lambda$ = ',num2str(round(mean(PRMS(:,1)),2))],'Interpreter','latex')
ylabel('Replication Rate, $\lambda$','Interpreter','latex')
xlabel('Number of Iterations','Interpreter','latex')
nexttile
plot(PRMS(:,2),'LineWidth',2,'Color',[.75,0,0])
grid on
legend(['$\bar\rho$ = ',num2str(round(mean(PRMS(:,2)),2))],'Interpreter','latex')
ylabel('Transmission Bias, $\rho$','Interpreter','latex')
xlabel('Number of Iterations','Interpreter','latex')
T.Title.String = 'Simulated Data';
shg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%         Parameter Estimation: Simulated Data         %%%%%%%%%%
% Load dataset: Weak
load('propagon_data_raw/Weak.mat','propagon_data','sampling_times')

% Setup initial parameter values to begin parameter estimation.
Prms0 = [0.1,0.5]; % Initial parameter estimates [lambda, rho]
Cn    = sqrt([0.1406,0;0,0.0156]); % Initial covariance matrix
num_itrs = 6*10^3;

% Perform parameter estimation using Adaptive Metropolis Algorithm
[PRMS,~] = am_alg(propagon_data,sampling_times,Prms0,Cn,num_itrs);

% Plot adaptive Metropolis chain
figure
L = tiledlayout(2,1);
nexttile
plot(PRMS(:,1),'LineWidth',2,'Color',[0,0,.75])
grid on
legend(['$\bar\lambda$ = ',num2str(round(mean(PRMS(:,1)),2))],'Interpreter','latex')
ylabel('Replication Rate, $\lambda$','Interpreter','latex')
xlabel('Number of Iterations','Interpreter','latex')
nexttile
plot(PRMS(:,2),'LineWidth',2,'Color',[.75,0,0])
grid on
legend(['$\bar\rho$ = ',num2str(round(mean(PRMS(:,2)),2))],'Interpreter','latex')
ylabel('Transmission Bias, $\rho$','Interpreter','latex')
xlabel('Number of Iterations','Interpreter','latex')
L.Title.String = 'Experimental Data';

set(gca,'TickLabelInterpreter','latex','FontSize',16)
shg