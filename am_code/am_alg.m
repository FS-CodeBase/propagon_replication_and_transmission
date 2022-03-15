function [PRMS,LL,Cmat] = am_alg(propagon_data,sampling_times,Prms0,Cn,nItrs)
% FUNCTION: AM_ALG 
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
% 
% DESCRIPTION: Applies the Adaptive Metropolis Algorithm to fit the asymmetric 
% % transmission of propagons model to the input propagon_data.
% 
% INPUTS:
% % propagon_data: aggregate data in a cell array.
% % sampling_times: hold the times that each aggreagate 
% % Prms0: initial parameter values
% % Cn: initial covariance matrix
% % nItrs: Number of iterations of the AM algorithm
% 
% OUTPUT:
% % PRMS: Parameter estimates by iterations
% % LL: Loglikelihood of parameter estimates
% % Cmat: Covariance matrix by iterations
%

% Load log version of aggreage dynamics model% [Yn,yn,cn,log_cn,log_Yn]
MU  = mean(propagon_data{1}); % Initial average distribution of aggregates
SIG = 1; % Set to model low variance of aggregates.
[~,~,cn,log_cn,log_Yn] = load_model_solutions(MU,SIG);

% Create anonymous function for Loglikelihood function of Data given
% aggregate replication dynamics (lambda).
LL =@(P) LogL(P, sampling_times, log_Yn, log_cn, cn, propagon_data);

% Pre-allocate space for lambda and rho estimates
% PRMS = [lambda rho];
PRMS = [Prms0;zeros(nItrs,2)];

PrmsPre = Prms0;
LLpre = LL(PrmsPre);

sD = 1; eD = 10^-6;
Id = eye(2).*(diag(Cn)>0);

dnon_adapt = 1501;
Cmat = zeros(2,2,nItrs-dnon_adapt);
Q = ones(2);
Q((~diag(Cn))==1,:)=0;
Q(:,(~diag(Cn))==1)=0;

% Adaptive Metropolis Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NON-ADAPTIVE STEP
for itr = 2:(dnon_adapt)
    PrmsTmp = abs(mvnrnd(PrmsPre,Cn,1));
    PrmsTmp([false PrmsTmp(2)>0.5]) = abs(1-PrmsTmp([false PrmsTmp(2)>0.5]));
    LLtmp = LL(PrmsTmp);
    lhs = log(rand())+LLpre;
    rhs = LLtmp;
    if lhs < rhs
        PRMS(itr,:) = PrmsTmp; 
        LLpre   = LLtmp;
        PrmsPre = PrmsTmp;
    else
        PRMS(itr,:) = PrmsPre;  
    end
end

% Adaptive Step
Cn = cov(PRMS(1:dnon_adapt,:));    % Covariance up to this point
muPrmsPre = mean(PRMS(1:dnon_adapt,:)); % Mean parameter values up to this point
for itr = (dnon_adapt+1):(nItrs+1)
    PrmsTmp = abs(mvnrnd(PrmsPre,Cn+eD*Id,1));
    
    % If rho>0.5, reflect about 0.50 so that rho in [0,0.5].
    PrmsTmp([false PrmsTmp(2)>0.5]) = 1-PrmsTmp([false PrmsTmp(2)>0.5]);
    
    LLtmp = LL(PrmsTmp);
    lhs = log(rand())+LLpre;
    rhs = LLtmp;
    if lhs < rhs
        PRMS(itr,:) = PrmsTmp; 
        LLpre   = LLtmp;
        PrmsPre = PrmsTmp;
    else
        PRMS(itr,:) = PrmsPre;
    end

    % Update Covariance Matrix using new parameter stamples
    muPrmsCur = (itr)/(itr+1)*muPrmsPre+1/(itr+1)*PRMS(itr,:);

    Cn = covar(Cn,itr,muPrmsPre,muPrmsCur,PRMS(itr,:),sD);
    Cn = Cn.*Q;
    Cmat(:,:,itr-dnon_adapt) = Cn;
    muPrmsPre = muPrmsCur;
end