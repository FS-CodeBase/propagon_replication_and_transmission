function Cnp1 = covar(Cn,n,barXp,barXn,Xn,sd)
% FUNCTION/SCRIPT: COVAR
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
% 
% DESCRIPTION
% % Updates the covariance matrix based on a new data point
% 
% INPUTS 
% % Cn: Current covariance matrix
% % n: number of samples to compute Cn
% % barXp: Previous mean parameter values n-1 samples
% % barXn: Current mean parameter values  n samples
% % Xn: Current parameter estimates
% % sd: Design paramter values
% 
% OUTPUT
% % Cnp1: Updated covariance matrix 
%
Cnp1 = (n-1)/n*Cn+sd/n*(n*((barXp')*barXp)-(n+1)*((barXn')*barXn)+((Xn')*Xn));