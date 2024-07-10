This folder contains the Matlab scripts and functions and OpenBUGS scripts, as well the necessary
input files, to implement the parameter estimation and uncertainty and sensitity analyses described
in Gubbins "Using the basic reproduction number to quantify transmission and identify data gaps
for epizootic haemorrhagic disease, an emerging vector-borne disease in Europe"

-------------------
MATLAB REQUIREMENTS
-------------------

The Matlab scripts/functions were run using Matlab version 2020b. The parameter estimation methods
require the Statistics and Machine Learning and Parallel Computing toolboxes. However, they can be
easily adapted to run without the Parallel Computing toolbox by changing any "parfor" loops to "for"
loops in the functions with names beginning "ParEst". Other scripts/functions do not require any
toolboxes.


-------------------------------------------------
MATLAB SCRIPTS/FUNCTIONS FOR PARAMETER ESTIMATION
-------------------------------------------------

estimateViraemiaInCattle.m - implements the adaptive Metropolis scheme to estimate the mean and shape
                             parameter for the duration of viraemia in cattle assuming no strain
                             dependence; the data are included in this file
estimateViraemiaInCattleToo.m - implements the adaptive Metropolis scheme to estimate the mean and shape
                                parameter for the duration of viraemia in cattle for different models
                                of their dependence on EHDV strain; the data are included in this file
estimateViraemiaDeer.m - implements the adaptive Metropolis scheme to estimate the mean and shape
                         parameter for the duration of viraemia in deer assuming no strain dependence;
                         the data are included in this file
Lhood_Viraemia.m - computes the log likelihood and prior for the input parameters for the duration of
                   viraemia assuming no strain dependence
LhoodToo_Viraemia.m - computes the log likelihood and prior for the input parameters for the duration of
                      viraemia for models making different assumptions about strain dependence of the
                      parameters

ParEst_EIP.m - loads the data and implements the adaptive Metropolis scheme to estimate EIP parameters for
               different models of strain dependence
Lhood_EIP.m - computes the log likelihood and prior for the input parameters for the EIP for the different
              models


-----------------------------------------
OpenBUGS SCRIPTS FOR PARAMETER ESTIMATION
-----------------------------------------

The following scripts estimate the parameters described in the file name using OpenBUGS. Each file
includes the data in the format used by the package.

- estimateBitingRate.txt
- estimateMortalityRate.txt
- estimateCaseFatalityInDeer.txt
- estimatePrH2V.txt


-------------------------------------------------------
MATLAB SCRIPTS FOR UNCERTAINTY AND SENSITIVITY ANALYSES
-------------------------------------------------------

samplePosteriors.m - loads the MCMC samples and generates parameter sets drawn from the joint posterior
                     distributions

computeR0VsTemperature.m - calculates R0 at different temperatures for parameter values sampled from
                           their joint posterior distributions; this calls samplePosteriors.m

computeSobolIndices.m - computes the first-order and total Sobol sensitivity indices for R0 at different
                        temperatures based on parameters sampled from their joint posterior distributions;
                        this calls samplePosteriors.m and computeR0VsTemperature.m


----------
DATA FILES
----------

EIPForEHDVData.txt - data from experiments measuring the EIP in Culicoides sonorensis; columns are:
                     expt - 1-4, data from Tables 1-4 in Ruder et al. (2015) J Med Entomol 52(5) 1050-1059
   	                    5, data from Figure 1 in Wittmann et al. (2002) Med Vet Entomol 16, 147-156
                     route - route of infection: 1-membrane; 2-viraemic deer
                     sero - EHDV serotype (1, 2 or 7)
                     T - rearing temperature (oC)
                     dpf - days post feeding
                     n - no. midges tested
                     nPH - no. midges with a fully disseminated infection

The following files contain the MCMC samples for the parameters described in the file name:
- CulicoidesParameters.txt (cols: 1-2 mortality rate (mu0, mu1); 3-4 biting rate (a0, T0)
- DurationOfViraemiaInCattle.txt (cols: 1 shape; 2 mean)
- DurationOfViraemiaInDeer.txt (cols: 1 shape; 2 mean)
- EIPParameters.txt (cols: 1-4 threshold temperature for each strain; 5-8 replication rate for each
                     strain; 9 shape; note: these are only the parameters used to compute R0, not the
                     full posterior samples from parameter estimation)
