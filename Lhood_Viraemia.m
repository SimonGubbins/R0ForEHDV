function [logL, prior]=Lhood_Viraemia(par,D)
%
% [logL, prior]=Lhood_Viraemia(par,D)
%
% Matlab function to compute the log likelihood and prior for the duration
% of viraemia in cattle or deer infected with EHDV
%
% Inputs:
% par - vector of parameters
% D - data on duration of viraemia; columns are:
%     1 - serotype
%     2 - minimum duration (days)
%     3 - maximum duration (days)
%     4 - no. cattle
%
% Outputs:
% logL - log likelihood for the parameters
% prior - log prior probability for the parameters

% Extract the parameters
sV=par(1);
muV=par(2);

% Compute the prior
prior=sum(log(exppdf(par,100)));

% Compute the log likelihood
logL=sum(D(:,4).*log((gamcdf(D(:,3),sV,muV./sV)-...
                      gamcdf(D(:,2),sV,muV./sV))));
