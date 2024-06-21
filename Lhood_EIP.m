function [logL, prior]=Lhood_EIP(par,D,mFlag)
%
% [logL, prior]=Lhood_EIP(par,D,mFlag)
%
% Matlab function for computing the log-likelihood for the temperature-
% dependent EIP model for EHDV.
%
% The EIP is assumed to follow a Gamma distribution with mean 1/v and
% variance 1/kv2, where v=v(T)=a*(T-Tmin) and k is the scale parameter.
% The analysis also allows for only a proportion (phi) of the vectors being
% competent to transmit virus.
%
% Inputs:
% par - vector containing the model parameters
% D - array containing the experimental data; columns are:
%     1 - experiment
%     2 - route of infection
%     3 - serotype
%     4 - temperature
%     5 - days post feeding when tested
%     6 - no. midges tested
%     7 - no. midges with a disseminated infection
% mFlag - flag indicating model to fit
%         1 - phi, Tmin, a common to all experiments
%         2 - Tmin, a common to all experiments, phi varies
%         3 - phi, a common to all experiments, Tmin varies
%         4 - phi, Tmin common to all experiments, a varies
%         5 - a common to all experiments, phi, Tmin vary
%         6 - Tmin common to all experiments, phi, a vary
%         7 - phi common to all experiments, Tmin, a vary
%         8 - phi, Tmin, a vary amongst experiments
%         9 - phi, Tmin, a vary amongst experiments (but Tmin and a for
%             experiments 3 & 4, which differ in route of infection, are
%             the same)
%
% Outputs:
% logL - log-likelihood for the parameter set
% prior - log prior

%==========================================================================
% EXTRACT THE PARAMETERS
% Set the number of experiments
nE=max(D(:,1));

% Extract and back-transform the parameters
if mFlag==1
    phi=repmat(par(1),nE,1);
    Tmin=repmat(exp(par(2)),nE,1);
    a=repmat(exp(par(3)),nE,1);
elseif mFlag==2
    phi=par(1:nE);
    Tmin=repmat(exp(par(nE+1)),nE,1);
    a=repmat(exp(par(nE+2)),nE,1);
    aP=exp(par(nE+3));
    bP=exp(par(nE+4));
elseif mFlag==3
    phi=repmat(par(1),nE,1);
    Tmin=exp(par(2:nE+1));
    a=repmat(exp(par(nE+2)),nE,1);
    aT=exp(par(nE+3));
    bT=exp(par(nE+4));
elseif mFlag==4
    phi=repmat(par(1),nE,1);
    Tmin=repmat(exp(par(2)),nE,1);
    a=repmat(exp(par(3:nE+2)),nE,1);
    aA=exp(par(nE+3));
    bA=exp(par(nE+4));
elseif mFlag==5
    phi=par(1:nE);
    Tmin=exp(par(nE+1:2*nE));
    a=repmat(exp(par(2*nE+1)),nE,1);
    aP=exp(par(2*nE+2));
    bP=exp(par(2*nE+3));
    aT=exp(par(2*nE+4));
    bT=exp(par(2*nE+5));
elseif mFlag==6
    phi=par(1:nE);
    Tmin=repmat(exp(par(nE+1)),nE,1);
    a=exp(par(nE+2:2*nE+1));
    aP=exp(par(2*nE+2));
    bP=exp(par(2*nE+3));
    aA=exp(par(2*nE+4));
    bA=exp(par(2*nE+5));
elseif mFlag==7
    phi=repmat(par(1),nE,1);
    Tmin=exp(par(2:nE+1));
    a=exp(par(nE+2:2*nE+1));
    aT=exp(par(2*nE+2));
    bT=exp(par(2*nE+3));
    aA=exp(par(2*nE+4));
    bA=exp(par(2*nE+5));
elseif mFlag==8
    phi=par(1:nE);
    Tmin=exp(par(nE+1:2*nE));
    a=exp(par(2*nE+1:3*nE));
    aP=exp(par(3*nE+1));
    bP=exp(par(3*nE+2));
    aT=exp(par(3*nE+3));
    bT=exp(par(3*nE+4));
    aA=exp(par(3*nE+5));
    bA=exp(par(3*nE+6));
elseif mFlag==9
    phi=par(1:nE);
    Tmin=exp(par(nE+1:2*nE-1));
    Tmin=[Tmin(1:2); Tmin(3); Tmin(3); Tmin(4)];
    a=exp(par(2*nE:3*nE-2));
    a=[a(1:2); a(3); a(3); a(4)];
    aP=exp(par(3*nE-1));
    bP=exp(par(3*nE));
    aT=exp(par(3*nE+1));
    bT=exp(par(3*nE+2));
    aA=exp(par(3*nE+3));
    bA=exp(par(3*nE+4));
end

% Extract the shape parameter
k=exp(par(end));
%==========================================================================

%==========================================================================
% COMPUTE THE PRIOR
% Compute the joint prior for the phi, Tmin and a
if mFlag==1
    prior=log(betapdf(phi(1),1,1))+...
          log(exppdf(Tmin(1),100))+...
          log(exppdf(a(1),100));
elseif mFlag==2
    prior=sum(log(betapdf(phi,aP,bP)))+...
          log(exppdf(Tmin(1),100))+...
          log(exppdf(a(1),100))+...
          log(exppdf(aP,100))+...
          log(exppdf(bP,100));
elseif mFlag==3
    prior=log(betapdf(phi(1),1,1))+...
          sum(log(gampdf(Tmin,aT,bT)))+...
          log(exppdf(a(1),100))+...
          log(exppdf(aT,100))+...
          log(exppdf(bT,100));
elseif mFlag==4
    prior=log(betapdf(phi(1),1,1))+...
          log(exppdf(Tmin(1),100))+...
          sum(log(gampdf(a,aA,bA)))+...
          log(exppdf(aA,100))+...
          log(exppdf(bA,100));
elseif mFlag==5
    prior=sum(log(betapdf(phi,aP,bP)))+...
          sum(log(gampdf(Tmin,aT,bT)))+...
          log(exppdf(a(1),100))+...
          log(exppdf(aP,100))+...
          log(exppdf(bP,100))+...
          log(exppdf(aT,100))+...
          log(exppdf(bT,100));
elseif mFlag==6
    prior=sum(log(betapdf(phi,aP,bP)))+...
          log(exppdf(Tmin(1),100))+...
          sum(log(gampdf(a,aA,bA)))+...
          log(exppdf(aP,100))+...
          log(exppdf(bP,100))+...
          log(exppdf(aA,100))+...
          log(exppdf(bA,100));
elseif mFlag==7
    prior=log(betapdf(phi(1),1,1))+...
          sum(log(gampdf(Tmin,aT,bT)))+...
          sum(log(gampdf(a,aA,bA)))+...
          log(exppdf(aT,100))+...
          log(exppdf(bT,100))+...
          log(exppdf(aA,100))+...
          log(exppdf(bA,100));
elseif mFlag==8
    prior=sum(log(betapdf(phi,aP,bP)))+...
          sum(log(gampdf(Tmin,aT,bT)))+...
          sum(log(gampdf(a,aA,bA)))+...
          log(exppdf(aP,100))+...
          log(exppdf(bP,100))+...
          log(exppdf(aT,100))+...
          log(exppdf(bT,100))+...
          log(exppdf(aA,100))+...
          log(exppdf(bA,100));
elseif mFlag==9
    prior=sum(log(betapdf(phi,aP,bP)))+...
          sum(log(gampdf(Tmin([1:3 5]),aT,bT)))+...
          sum(log(gampdf(a([1:3 5]),aA,bA)))+...
          log(exppdf(aP,100))+...
          log(exppdf(bP,100))+...
          log(exppdf(aT,100))+...
          log(exppdf(bT,100))+...
          log(exppdf(aA,100))+...
          log(exppdf(bA,100));
end

% Add the prior for the shape parameter
prior=prior+log(exppdf(k,100));
%==========================================================================

%==========================================================================
% COMPUTE THE LOG LIKELIHOOD
% Calculate the reciprocal of the mean EIP for the sample temperatures
nu=max(0,a(D(:,1)).*(D(:,4)-Tmin(D(:,1))));

% Compute the probability that a sampled midge is positive
prPos=phi(D(:,1)).*gamcdf(D(:,5),k,1./(k.*nu));

% Compute the log-likelihood
logL=sum(log(binopdf(D(:,7),D(:,6),prPos)));
%==========================================================================
