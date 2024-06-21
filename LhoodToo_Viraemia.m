function [logL, prior]=LhoodToo_Viraemia(par,D,mFlag)
%
% [logL, prior]=LhoodToo_Viraemia(par,D,mFlag)
%
% Matlab function to compute the log likelihood and prior for the duration
% of viraemia in cattle or deer infected with EHDV
%
% This version allows parameters to vary amongst strains
%
% Inputs:
% par - vector of parameters
% D - data on duration of viraemia; columns are:
%     1 - serotype
%     2 - minimum duration (days)
%     3 - maximum duration (days)
%     4 - no. cattle
% mFlag - flag indicating model to fit
%         1 - shape and mean common to all strains
%         2 - shape varies amongst strains and mean common to all strains
%         3 - shape common to all strains and mean varies amongst strains
%         4 - shape and mean vary amongst strains
%
% Outputs:
% logL - log likelihood for the parameters
% prior - log prior probability for the parameters

% Extract the parameters
if mFlag==1
    sV=par(1)*ones(2,1);
    muV=par(2)*ones(2,1);
elseif mFlag==2
    sV=par(1:2);
    muV=par(3)*ones(2,1);
elseif mFlag==3
    sV=par(1)*ones(2,1);
    muV=par(2:3);
elseif mFlag==4
    sV=par(1:2);
    muV=par(3:4);
end

% Compute the prior
prior=sum(log(exppdf(par,100)));

% Compute the log likelihood
logL=sum(D(:,4).*log((gamcdf(D(:,3),sV(D(:,1)),muV(D(:,1))./sV(D(:,1)))-...
                      gamcdf(D(:,2),sV(D(:,1)),muV(D(:,1))./sV(D(:,1))))));
