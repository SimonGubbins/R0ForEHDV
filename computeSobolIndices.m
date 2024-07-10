function [si,sTi]=computeSobolIndices(T,spp,sFlag,nSamp)
%
% [si,sTi]=computeSobolIndices(T,spp,sFlag,nSamp)
%
% Matlab function to compute first-order and total Sobol indices for R0
% at a given temperature for different strains of epizootic haemorrhagic
% disease virus in cattle and deer using random sampling from the (joint)
% posterior distributions
%
% This uses the Monte Carlo methods described in Azzani et al.
% Environmental Modelling and Software 144, 105167 (2021)
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
% Sample from the joint posterior distribution ...
parsA=samplePosteriors(spp,sFlag,nSamp);

% ... and compute the corresponding R0
R0A=computeR0VsTemperature(spp,parsA,T);

% Reample from the joint posterior distribution ...
parsB=samplePosteriors(spp,sFlag,nSamp);

% ... and compute the corresponding R0
R0B=computeR0VsTemperature(spp,parsB,T);
%==========================================================================

%==========================================================================
% COMPUTE THE SOBOL INDICES
% Create arrays to store the indices
si=NaN(size(parsA,2),length(T));
sTi=NaN(size(parsA,2),length(T));

% For each parameter ...
for j=1:size(parsA,2)

% Create paramter sets with all columns the same except current parameter
% which is set to that in the other set
    parsAB=parsA;
    parsAB(:,j)=parsB(:,j);
    parsBA=parsB;
    parsBA(:,j)=parsA(:,j);

% Compute R0 for the new parameter sets
    R0AB=computeR0VsTemperature(spp,parsAB,T);
    R0BA=computeR0VsTemperature(spp,parsBA,T);

% Compute the first-order Sobol index for the parameter
    si(j,:)=2.*sum((R0BA-R0B).*(R0A-R0AB),1)./....
               sum((R0A-R0B).^2+(R0BA-R0AB).^2,1);

% Compute the Sobol index for the total effect of the parameter
    sTi(j,:)=sum((R0B-R0BA).^2+(R0A-R0AB).^2,1)./...
             sum((R0A-R0B).^2+(R0BA-R0AB).^2,1);

end
%==========================================================================
