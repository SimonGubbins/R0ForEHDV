function [si,sTi]=computeSobolIndices(T,spp,sFlag,nSamp)
%
% [si,sTi]=computeSobolIndices_Cattle(T,sFlag,nSamp)
%
% Matlab function to compute first-order and total Sobol indices for R0
% at a given temperature for different strains of epizootic haemorrhagic
% disease virus in cattle and deer using random sampling from the (joint)
% posterior distributions
%
% Inputs:
% T - vector temperatures at which to compute R0 (and, hence, the indices)
% spp - string indication species being considered ('Cattle' or 'Deer')
% sFlag - flag indicating strain for which to compute R0:
%         1 - EHDV-1 (USA, deer)
%         2 - EHDV-2 (USA, deer)
%         3 - EHDV-7 (Israel, cattle)
%         4 - EHDV-1 (South Africa, unknown)
% nSamp - number of samples to draw
%
% Outputs:
% si - first-order Sobol indices for each parameter at each temperature
% sTi - Sobol indices for total effect of each parameter at each
%       temperature

%==========================================================================
% SAMPLE THE PARAMETERS AND COMPUTE R0
% Sample from the joint posterior distribution and compute the
% corresponding R0s
pars=samplePosteriors(spp,sFlag,nSamp);
R0=computeR0VsTemperature(spp,pars,T);

% Compute the total variance for R0
f0=sum(R0,1)./nSamp;
V=sum(R0.^2,1)./nSamp-(f0.^2);
%==========================================================================

%==========================================================================
% COMPUTE THE SOBOL INDICES
% Create arrays to store the indices
si=NaN(size(pars,2),length(T));
sTi=NaN(size(pars,2),length(T));

% Resample the parameters and compute the corresponding R0s
parsR=samplePosteriors(spp,sFlag,nSamp);
R0R=computeR0VsTemperature(spp,parsR,T);

% For each parameter ...
for j=1:size(pars,2)

% Create a paramter set with all columns except the current parameter set
% to the resampled values
    pars2=parsR;
    pars2(:,j)=pars(:,j);

% Compute R0 for the new parameter set
    R02=computeR0VsTemperature(spp,pars2,T);

% Compute the first-order Sobol index for the parameter
    Vi=sum(R0.*R02,1)./nSamp-(f0.^2);
    si(j,:)=Vi./V;

% Compute the Sobol index for the total effect of the parameter
    Vni=sum(R0R.*R02,1)./nSamp-(f0.^2);
    sTi(j,:)=1-Vni./V;

end
%==========================================================================
