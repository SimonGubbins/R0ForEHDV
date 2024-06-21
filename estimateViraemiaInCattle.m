function estimateViraemiaInCattle(seeds,nsamp,nburnin,nthin)
%
% estimateViraemiaInCattle(seeds,nsamp,nburnin,nthin)
%
% Matlab function to fit a gamma distribution to data on the duration of
% EHDV viraemia in cattle using Bayesian MCMC methods
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

% Data for duration of viraemia in cattle. Columns in each array are:
% 1 - serotype
% 2 - minimum duration (days)
% 3 - maximum duration (days)
% 4 - no. cattle
 
% 1) Data from Gard & Melville 1992 (copied from Table 6 in Savini et al.
% 2011)
D1=[2 0  7  46;
    2 7  14 24;
    2 14 21 3;
    2 21 28 4;
    2 28 35 1;
    5 0  7  31;
    5 7  14 3;
    5 21 28 1;
    5 28 35 1;
    7 7  14 1;
    8 0  7  12];

% 2) Data from Aradaib et al. 1994 (their Table 1)
D2=[1 14 18 1;
    2 7  11 1];

% 3) Data from Eschbaumer et al. 2012 (their Figure 4)
D3=[7 30 38 1;
    7 4 6 1;
    7 3 5 1;
    7 8 13 1];

% 4) Data from Ruder et al. 2015 (their Table III)
D4=[7 6 11 1;
    7 4 9  2];

% 5) Data from Spedicato et al. 2024 (their Table I)
D5=[8 1 7 1;
    8 9 19 1;
    8 14 Inf 1;
    8 4 12 1];

% Merge the data sets
D=[D1; D2; D3; D4; D5];

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
save('EHDV_DurationOfViraemiaInCattle_MCMCSamples',...
     'ParSamp','nburnin','nsamp','seeds')
