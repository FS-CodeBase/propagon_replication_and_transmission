function res = rem_outliers_iqr(in_counts)
% FUNCTION: rem_outliers_iqr
% % Author: Fabian Santiago
% % E-mail: FabianSantiago707@gmail.com
%
% DESCRIPTION 
% % Removes outliers using the inner quartile range (IQR). Counts that
% % fall outside of the range (Q1-1.5*IQR, Q3+1.5*IQR) are removed.
%

% Determine the first three quantiles
Q = quantile(in_counts,3);

% Compute IQR
IQR = Q(3)-Q(1);

% Remove outliers by keeping non-outliers
res = in_counts((in_counts>(Q(1)-1.5*IQR))&(in_counts<(Q(3)+1.5*IQR)));