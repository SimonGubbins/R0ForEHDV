function ParEst_EIP(mFlag,seeds,nsamp,nburnin,nthin)
%
% ParEst_EIP(mFlag,seeds,nsamp,nburnin,nthin)
%
% Matlab function to find implement a Bayesian MCMC scheme to estimate the
% extrinsic incubation period for EHDV in Culicoides sonorensis using data
% presented in Wittmann et al. 2002 Med Vet Entomol 16, 147-156 and Ruder
% et al. 2015 J Med Entomol 52, 1050-1059
%
% Inputs:
% mFlag - flag indicating model to fit:
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
% route - string indicating source of infection ('Membrane' or 'Deer')
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

% Load the experimental data. Columns are:
% 1 - experiment (Ruder-1, Wittmann-2)
% 2 - route of infection (membrane-1, deer-2)
% 3 - serotype
% 4 - temperature
% 5 - days post feeding when tested
% 6 - no. midges tested
% 7 - no. midges with a high titre (assumed indicative of a disseminated
%     infection)
D=readcell('EIPForEHDVData.txt','NumHeaderLines',1);
D=cell2mat(D);

% Set the number of parameters
nE=max(D(:,1));
if mFlag==1
    npar=3+1;
elseif mFlag==2 || mFlag==3 || mFlag==4
    npar=(nE+2)+2+1;
elseif mFlag==5 || mFlag==6 || mFlag==7
    npar=2*(nE+2)+1+1;
elseif mFlag==8
    npar=3*(nE+2)+1;
elseif mFlag==9
    npar=(nE+2)+2*(nE-1+2)+1;
end

% Create the arrays storing the output for the chain
ParSamp=cell(1,nchains);

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
        if mFlag==1
            phi=unifrnd(0.5,1);
            Tmin=unifrnd(15,20);
            a=unifrnd(0,0.01);
            par=[phi; log(Tmin); log(a)];
        elseif mFlag==2
            phi=unifrnd(0.5,1,nE,1);
            p=betafit(phi);
            aP=p(1);
            bP=p(2);
            Tmin=unifrnd(15,20);
            a=unifrnd(0,0.01);
            par=[phi; log(Tmin); log(a); log([aP; bP])];
        elseif mFlag==3
            phi=unifrnd(0.5,1);
            Tmin=unifrnd(15,20,nE,1);
            p=gamfit(Tmin);
            aT=p(1);
            bT=p(2);
            a=unifrnd(0,0.01);
            par=[phi; log(Tmin); log(a); log([aT; bT])];
        elseif mFlag==4
            phi=unifrnd(0.5,1);
            Tmin=unifrnd(15,20);
            a=unifrnd(0,0.01,nE,1);
            p=gamfit(a);
            aA=p(1);
            bA=p(2);
            par=[phi; log(Tmin); log(a); log([aA; bA])];
        elseif mFlag==5
            phi=unifrnd(0.5,1,nE,1);
            p=betafit(phi);
            aP=p(1);
            bP=p(2);
            Tmin=unifrnd(15,20,nE,1);
            p=gamfit(Tmin);
            aT=p(1);
            bT=p(2);
            a=unifrnd(0,0.01);
            par=[phi; log(Tmin); log(a); log([aP; bP; aT; bT])];
        elseif mFlag==6
            phi=unifrnd(0.5,1,nE,1);
            p=betafit(phi);
            aP=p(1);
            bP=p(2);
            Tmin=unifrnd(15,20);
            a=unifrnd(0,0.01,nE,1);
            p=gamfit(a);
            aA=p(1);
            bA=p(2);
            par=[phi; log(Tmin); log(a); log([aP; bP; aA; bA])];
        elseif mFlag==7
            phi=unifrnd(0.5,1);
            Tmin=unifrnd(15,20,nE,1);
            p=gamfit(Tmin);
            aT=p(1);
            bT=p(2);
            a=unifrnd(0,0.01,nE,1);
            p=gamfit(a);
            aA=p(1);
            bA=p(2);
            par=[phi; log(Tmin); log(a); log([aT; bT; aA; bA])];
        elseif mFlag==8
            phi=unifrnd(0.5,1,nE,1);
            p=betafit(phi);
            aP=p(1);
            bP=p(2);
            Tmin=unifrnd(15,20,nE,1);
            p=gamfit(Tmin);
            aT=p(1);
            bT=p(2);
            a=unifrnd(0,0.01,nE,1);
            p=gamfit(a);
            aA=p(1);
            bA=p(2);
            par=[phi; log(Tmin); log(a); log([aP; bP; aT; bT; aA; bA])];
        elseif mFlag==9
            phi=unifrnd(0.5,1,nE,1);
            p=betafit(phi);
            aP=p(1);
            bP=p(2);
            Tmin=unifrnd(15,20,nE-1,1);
            p=gamfit(Tmin);
            aT=p(1);
            bT=p(2);
            a=unifrnd(0,0.01,nE-1,1);
            p=gamfit(a);
            aA=p(1);
            bA=p(2);
            par=[phi; log(Tmin); log(a); log([aP; bP; aT; bT; aA; bA])];
        end
        
% Generate a shape parameter
        k=unifrnd(1,4);
        par=[par; log(k)];

% Compute the log-likelihood
        [CurrL, prior]=Lhood_EIP(par,D,mFlag);

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
        [NewL, prior_new]=Lhood_EIP(par_new,D,mFlag);

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

%==========================================================================
% COMPUTE DIC AND pD
% Compute the deviance for each sample
Dev=[];
PS=[];
for chain=1:nchains
    Dev=[Dev; -2*ParSamp{chain}(:,end)];
    PS=[PS; ParSamp{chain}(:,1:end-2)];
end

% Compute the mean deviance
Dbar=mean(Dev);

% Compute the deviance at the posterior mean for the parameters
Dhat=-2*Lhood_EIP(mean(PS,1)',D,mFlag);

% Compute the DIC
DIC=2*Dbar-Dhat;

% Compute the effective number of parameters
pD=Dbar-Dhat;
%==========================================================================

% Save the outputs
save(['EHDV_MCMCSamples_Model' num2str(mFlag)],...
     'ParSamp','nburnin','nsamp','seeds','DIC','pD')
