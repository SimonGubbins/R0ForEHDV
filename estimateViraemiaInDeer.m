function estimateViraemiaInDeer(seeds,nsamp,nburnin,nthin)
%
% estimateViraemiaInDeer(seeds,nsamp,nburnin,nthin)
%
% Matlab function to fit a gamma distribution to data on the duration of
% EHDV viraemia in white-tailed deer using Bayesian MCMC methods
%
% Inputs:
% seeds - vector of seeds to use for the random number generator for each
%         chain
% nsamp - number of samples to use when estimating parameters
% nburnin - number of samples to discard before estimating parameters
% nthin - number of samples by which to thin each chain
%
% Outputs: none (N.B. All chains are saved rather than provided as output
% arguments)

%==========================================================================
% PREPARATORY STUFF
% Set the number of chains
nchains=length(seeds);

% Data for duration of viraemia in white-tailed deer. Columns in each
% array are:
% 1 - serotype
% 2 - minimum duration (days)
% 3 - maximum duration (days)
% 4 - dead (1) or alive (0) at end of experiment
 
% 1) Data from Quist et al. 1997 (their Table 1)
D1=[2 6  Inf 1
    2 4  Inf 1
    2 11 Inf 1
    2 20 Inf 0
    2 20 Inf 0
    2 20 Inf 0
    2 4  8   0
    2 20 Inf 0
    2 6  Inf 1
    2 17 21  0
    2 10 Inf 1
    2 17 21  0
    2 13 17  0
    2 13 17  0
    2 20 Inf 0
    2 20 Inf 0];

% 2) Data from Ruder et al. 2012 (their Table 1)
D2=[7 2  Inf 1
    7 36 39  0
    7 4  Inf 1
    7 4  Inf 1
    7 4  Inf 1
    7 29 32  0];

% 3) Data from Ruder et al. 2016 (their Table I)
D3=[6 9  Inf 1
    6 9  Inf 1
    6 14 Inf 0
    6 16 Inf 0
    6 1  Inf 1];

% Merge the data sets
D=[D1(:,1:3); D2(:,1:3); D3(:,1:3)];
D=[D ones(size(D,1),1)];

% Set the number of parameters
npar=2;

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);
%==========================================================================

% For each chain ...
parfor chain=1:nchains

%==========================================================================
% INITIALISE THE CHAIN
% Initialise the random number generator
    rng(seeds(chain),'twister');

% Set the initial scaling factor for the proposal distribution
    sf=(2.38.^2)/npar;
    SIG=eye(npar);

% Set the counter for the number of accepted samples
    n_accept=0;

% Create the arrays storing the output for the chain
    ParSampC=zeros(nsamp/nthin,npar+2);
    iter=1;

% Generate the initial parameters for the chain, ensuring they generate a
% finite log likelihood and prior
    disp('Initialising chain')
    CurrL=NaN;
    prior=NaN;
    while ~isfinite(CurrL+prior)

% Generate an initial set of parameters
        sV=unifrnd(1,3);
        muV=unifrnd(5,15);
        par=[sV; muV];

% Compute the log-likelihood
        [CurrL, prior]=Lhood_Viraemia(par,D);

    end
%==========================================================================

%==========================================================================
% UPDATE THE PARAMETERS
% Sample parameter space
    disp('Sampling parameter space')
    for samp=1:nsamp+nburnin

% Indicate what's going on
        disp(['Chain: ' num2str(chain) ', Sample: ' num2str(samp) ';'...
              ' Accept: ' num2str(100*n_accept/samp,3) '%'])

% Update the variance-covariance matrix for the proposal distribution
        if samp<=nburnin && (samp<=2*npar || n_accept==0)
            SIGp=0.01*eye(npar);
        else
            SIGp=sf.*(SIG+0.01*eye(npar));
        end

% Generate the new set of probabilities
        par_new=par+mvnrnd(zeros(1,length(par)),SIGp)';

% Compute the log likelihood and prior for the new parameter set
        [NewL, prior_new]=Lhood_Viraemia(par_new,D);

% Test whether to accept the new parameter set
        u=unifrnd(0,1);
        if isfinite(NewL+prior_new) && ...
           u<min(1,exp((NewL+prior_new)-(CurrL+prior)))

% Update the counter
            n_accept=n_accept+1;

% Update the covariance matrix for the proposal distribution
            if n_accept==1
                pbar=mean([par par_new],2);
                SIG=cov([par'; par_new']);
            elseif samp<=nburnin && n_accept>1
                pbar_prev=pbar;
                pbar=(n_accept./(n_accept+1)).*pbar_prev+...
                     (1./(n_accept+1)).*par_new;
                SIG=((n_accept-1)./n_accept).*SIG+...
                    (1./n_accept).*(n_accept.*(pbar_prev*pbar_prev')-...
                                    (n_accept+1).*(pbar*pbar')+...
                                    (par_new*par_new'));
            end

% Update the chain
            CurrL=NewL;
            prior=prior_new;
            par=par_new;

        end

% Every one hundred samples during burn-in, tune the scaling factor
% for the proposal distribution to ensure an acceptance rate of 20-40%
        if samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp<0.2
            sf=sf/2;
        elseif samp<=nburnin && mod(samp+1,100)==1 && n_accept/samp>0.4
            sf=2*sf;
        end
%==========================================================================

%==========================================================================
% STORE THE OUTPUT
% After burn in, save iterations of the chain, thinning as specified
        if nthin==1
            ParSampC(samp,:)=[par' prior CurrL];
        elseif samp>nburnin && mod(samp,nthin)==1
            ParSampC(iter,:)=[par' prior CurrL];
            iter=iter+1;
        end
%==========================================================================

    end

% Store the chain
    ParSamp{chain}=ParSampC;

end

% Save the outputs
save('EHDV_DurationOfViraemiaInDeer_MCMCSamples',...
     'ParSamp','nburnin','nsamp','seeds')
