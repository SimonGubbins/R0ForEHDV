function [si,sTi]=computeSobolIndices(T,spp,sFlag,nSamp)
%
% [si,sTi]=computeSobolIndices(T,spp,sFlag,nSamp)
%
% Matlab function to compute first-order and total Sobol indices for R0
% at a given temperature for different strains of epizootic haemorrhagic
% disease virus in cattle or deer using random sampling from the (joint)
% posterior distributions
%
% This uses the Monte Carlo methods described in Azzani et al.
% Environmental Modelling and Software 144, 105167 (2021) and treats
% jointly sampled parameters as groups rather than independently when
% resampling
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

% Specify the parameter groupings. Vector elements are:
% 1 - probability of transmission (vector to host)
% 2 - probability of transmission (host to vector)
% 3 - vector to host ratio
% 4 - biting rate parameters (a0, T0)
% 5 - duration of viraemia (mean, shape) and disease associated mortality
% 6 - EIP parameters (replication rate, threshold temperature, shape)
% 7 - mortality rate parameters (mu0, mu1)
if strcmp(spp,'Cattle')
    parGrp={1, 2, 3, 4:5, 6, 7:9, 10:11};
elseif strcmp(spp,'Deer')
    parGrp={1, 2, 3, 4:5, 6:8, 9:11, 12:13};
end

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
si=NaN(length(parGrp),length(T));
sTi=NaN(length(parGrp),length(T));

% For each parameter group...
for j=1:length(parGrp)

% Create paramter sets with all columns the same except current parameter
% which is set to that in the other set
    parsAB=parsA;
    parsAB(:,parGrp{j})=parsB(:,parGrp{j});
    parsBA=parsB;
    parsBA(:,parGrp{j})=parsA(:,parGrp{j});

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
