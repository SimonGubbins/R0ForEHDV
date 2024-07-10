function pars=samplePosteriors(spp,sFlag,nSamp)
%
% pars=samplePosteriors(spp,sFlag,nSamp)
%
% Matlab function to generate replicated samples from the (joint) posterior
% distribution(s) for parameters used to compute the basic reproduction
% number for epizootic haemorrhagic disease virus
%
% Inputs:
% spp - string indication species being considered ('Cattle' or 'Deer')
% sFlag - flag indicating strain for which to compute R0:
%         1 - EHDV-1 (USA, deer)
%         2 - EHDV-2 (USA, deer)
%         3 - EHDV-7 (Israel, cattle)
%         4 - EHDV-1 (South Africa, unknown)
% nSamp - number of samples to draw
%
% Outputs:
% pars - array of replicated, sampled parameters (rows are samples; see
%        below for columns)

% NOTE
% Columns in the pars array for cattle are:
% 1 - probability of transmission (vector to host)
% 2 - probability of transmission (host to vector)
% 3 - vector to host ratio
% 4 - biting rate (threshold temperature)
% 5 - biting rate (threshold temperature)
% 6 - mean duration of viraemia
% 7 - virus replication rate
% 8 - threshold temperature for replication
% 9 - shape parameter for EIP
% 10 - mortality rate parameter (mu0)
% 11 - mortality rate parameter (mu1)
%
% Columns in the pars array for deer are:
% 1 - probability of transmission (vector to host)
% 2 - probability of transmission (host to vector)
% 3 - vector to host ratio
% 4 - biting rate (threshold temperature)
% 5 - biting rate (threshold temperature)
% 6 - mean duration of viraemia
% 7 - shape parameter for duration of viraemia
% 8 - disease associated mortality rate
% 9 - virus replication rate
% 10 - threshold temperature for replication
% 11 - shape parameter for EIP
% 12 - mortality rate parameter (mu0)
% 13 - mortality rate parameter (mu1)

%==========================================================================
% CREATE THE POSTERIOR DISTRIBUTIONS
% Specify the posterior distributions to use for each parameter:
% 1 - probability of transmission (vector to host)
% 2 - probability of transmission (host to vector)
% 3 - vector to host ratio
% 4 - joint distribution for biting rate (slope and threshold temperature)
% 5 - joint distribution for mean duration of viraemia, shape parameter and
%     disease-associated mortality rate
% 6 - joint distribution for virus replication rate, threshold temperature
%     for replication and shape parameter for EIP
% 7 - joint distribution for mortality rate (mu0 and mu1)
pDist={'Beta',[7.38 2.12];
       'Beta',[1.36 47.8];
       'Gamma',[1.62 1774/1.62];
       [],[];
       [],[];
       [],[];
       [],[]};

% Load the MCMC samples for the duration of viraemia parameters and compute
% estimates for the disease-associated mortality rate
PS=load(['DurationOfViraemiaIn' spp '.txt']);
if strcmp(spp,'Cattle')
    pDist{5,1}=PS(:,2);
elseif strcmp(spp,'Deer')
    muD=PS(:,2);
    sD=PS(:,1);
    fD=betarnd(71.0,87.1,length(muD),1);
    dD=(sD./muD).*(1./((1-fD).^(1./sD))-1);
    pDist{5,1}=[muD sD dD];
end

% Load the MCMC samples for the Culicoides biting and mortality rates
PS=load('CulicoidesParameters.txt');
pDist{4,1}=PS(:,3:4);
pDist{7,1}=PS(:,1:2);

% Load the MCMC samples for the EIP parameters
PS=load('EIPParameters.txt');
pDist{6,1}=PS(:,[4+sFlag sFlag end]);
%==========================================================================

%==========================================================================
% SAMPLE FROM THE JOINT POSTERIOR DISTRIBUTION
% Create the array to store the sampled parameters
if strcmp(spp,'Cattle')==1
    pars=NaN(nSamp,11);
elseif strcmp(spp,'Deer')==1
    pars=NaN(nSamp,13);
end

% Sample the parameters:
% 1 - probability of transmission (vector to host)
% 2 - probability of transmission (host to vector)
% 3 - vector to host ratio
for j=1:3
    pars(:,j)=random(pDist{j,1},pDist{j,2}(1),pDist{j,2}(2),nSamp,1);
end

% Sample from the joint distribution for biting rate (slope and threshold
% temperature)
u=unidrnd(size(pDist{4,1},1),nSamp,1);
pars(:,4:5)=pDist{4,1}(u,:);

% Sample the mean duration of viraemia
u=unidrnd(size(pDist{5,1},1),nSamp,1);
if strcmp(spp,'Cattle')==1
    pars(:,6)=pDist{5,1}(u,:);
elseif strcmp(spp,'Deer')==1
    pars(:,6:8)=pDist{5,1}(u,:);
end

% Sample from the joint distribution for virus replication rate, threshold
% temperature for replication and shape parameter for EIP
u=unidrnd(size(pDist{6,1},1),nSamp,1);
if strcmp(spp,'Cattle')==1
    pars(:,7:9)=pDist{6,1}(u,:);
elseif strcmp(spp,'Deer')==1
    pars(:,9:11)=pDist{6,1}(u,:);
end

% Sample from the joint distribution for mortality rate (mu0 and mu1)
u=unidrnd(size(pDist{7,1},1),nSamp,1);
if strcmp(spp,'Cattle')==1
    pars(:,10:11)=pDist{7,1}(u,:);
elseif strcmp(spp,'Deer')==1
    pars(:,12:13)=pDist{7,1}(u,:);
end
%==========================================================================
