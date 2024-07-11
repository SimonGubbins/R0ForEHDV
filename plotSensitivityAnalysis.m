function plotSensitivityAnalysis(T,nSamp,nReps)
%
% plotSensitivityAnalysis(T,nSamp,nReps)
%
% Matlab function to plot sensitivity analysis of R0 as a function of
% temperature for different strains of epizootic haemorrhagic disease
% virus in cattle and deer. Specifically, it uses Monte Carlo methods to
% compute the first and total Sobol sensitivity indices.
%
% Inputs:
% T - vector of temperatures at which to compute the sensitivity measures
% nSamp - number of samples to draw from the joint posterion distribution
%         for the model parameters
% nReps - number of replicates for which to compute the sensitivity
%         indices
%
% Outputs: (none)

%==========================================================================
% PLOTTING STUFF
% Specify the parameter labels
parLab={'{\it{b}}', '\beta', '{\it{m}}', 'biting', 'host pars', 'EIP', ...
        'vec. mort.'};

% List the species
sppList={'Cattle'; 'Deer'};
sppLab={'cattle'; 'deer'};

% List the strains
sLab={'EHDV-1 (USA)',...
      'EHDV-2 (USA)',...
      'EHDV-7 (Israel)',...
      'EHDV-1 (unknown)'};

% Set the axis locations
wd=0.18;
ht=0.20;
ax_xpos=0.11:0.20:0.71;
ax_ypos=0.77:-0.23:0.08;

% Create a colour map
C=flip(colormap(pink),1);
C=[0.7*ones(128,3); C(2:2:256,:)];
%==========================================================================

% For each species ...
for spp=1:length(sppList)

% For each strain ...
    for s=1:length(sLab)
        disp([' ' sppLab{spp} ', ' sLab{s}])
    
%==========================================================================
% COMPUTE THE SOBOL INDICES
% Create arrays to store the Sobol indices
        s1=NaN(length(parLab),length(T),nReps);
        sT=NaN(length(parLab),length(T),nReps);

% For each replicate ...
        for r=1:nReps

% Compute the Sobol indices
            [s1(:,:,r),sT(:,:,r)]=computeSobolIndices(T,sppList{spp},...
                                                      s,nSamp);
        end

% Store the Sobol indices for each parameter
        s1=median(s1,3);
        sT=median(sT,3);

% Adjust any sensitivity measures that are NaNs to so they will be grey for
% plotting
        s1(isnan(s1))=-1;
        sT(isnan(sT))=-1;
%==========================================================================

%==========================================================================
% PLOT SENSITIVITY INDICES VS TEMPERATURE
% For each index ...
        for i=1:2

% Select the one to plot
            if i==1
                sensInd=s1;
            elseif i==2
                sensInd=sT;
            end

% Select the figure
            loc=2*(spp-1)+i;
            subplot('position',[ax_xpos(s) ax_ypos(loc) wd ht])
        
% Plot the sensitivity measures (note: ranges set so non-computed values
% appear as grey)
            imagesc(sensInd,[-1 1]);
            colormap(C)

% Add some dividing lines
            for j=1:length(parLab)-1
                line([0.5 length(T)+0.5],(j+0.5)*[1 1],'linestyle','-',...
                     'color','k','linewidth',0.25)
            end

% Add a colour bar
            if s==4
                colorbar;
                ax=colorbar;
                ax.FontSize=6;
                ax.YLim=[0 1];
                ax.YTick=0:0.2:1;
                ax.YTickLabel={'0','0.2','0.4','0.6','0.8','1.0'};
                ax.Label.FontSize=6;
                if i==1
                    ax.Label.String='first-order sensitivity index';
                elseif i==2
                    ax.Label.String='total sensitivity index';
                end
                ax.Position=[ax_xpos(s)+wd+0.01 ax_ypos(loc) 0.02 ht];
            end

% Make the axes look pretty
            ax=gca;
            ax.FontSize=6;
            ax.Box='on';
            ax.XTick=1:5:length(T);
            if loc==4
                ax.XTickLabel=T(1:5:length(T));
            else
                ax.XTickLabel={''};
            end
            ax.YTick=1:length(parLab);
            if s==1
                ax.YTickLabel=parLab;
            else
                ax.YTickLabel={''};
            end

% Label the x-axis
            if loc==4
                xlabel('temperature ({\circ}C)')
            end

% Indicate the species
            if s==1
                ylabel(sppLab{spp})
            end

% Indicate the strain
            if spp==1 && i==1
                title(sLab{s},'fontsize',6,'fontweight','normal')
            end
            
        end
%==========================================================================

    end
end

% Print to a file
print('-dtiff','-r300','..\SobolIndices_R0VsTemperatureForEHDV.tif')

% Tidy up
close('all')
