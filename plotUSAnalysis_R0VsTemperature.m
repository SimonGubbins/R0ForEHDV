function plotUSAnalysis_R0VsTemperature(T,sensMeas,nSamp,nReps)
%
% plotUSAnalysis_R0VsTemperature(T,sensMeas,nSamp,nReps)
%
% Matlab function to implement sensitivity analysis of R0 as a function of
% temperature for different strains of epizootic haemorrhagic disease
% virus in cattle and deer
%
% Inputs:
% T - vector of temperatures at which to compute the sensitivity measures
% sensMeas - string identifying sensitivity measure to compute: 'Sobol1'
%            for first-order Sobol indices; or 'SobolT' for total Sobol
%            indices
% nSamp - number of samples to draw from the joint posterion distribution
%         for the model parameters
% nReps - number of replicates for which to compute the sensitivity
%         measures
%
% Outputs: (none)

%==========================================================================
% PLOTTING STUFF
% Specify the parameter labels
parLab={{'{\it{b}}', '\beta', '{\it{m}}', '{\it{a}}_0', '{\it{T}}_0',...
         '1/{\it{r}}', '{\alpha}', '{\it{T}}_{min}', '\it{k}', ...
         '{\mu}_0', '{\mu}_1'}; ...
	    {'{\it{b}}', '\beta', '{\it{m}}', '{\it{a}}_0', '{\it{T}}_0',...
         '1/{\it{r}}', '{\it{n}}', '{\it{d}}', '{\alpha}', ...
         '{\it{T}}_{min}', '\it{k}', '{\mu}_0', '{\mu}_1'}};

% List the species
sppList={'cattle'; 'deer'};
sppListToo={'Cattle'; 'Deer'};

% List the strains
sLab={'EHDV-1 (USA)',...
      'EHDV-2 (USA)',...
      'EHDV-7 (Israel)',...
      'EHDV-1 (unknown)'};

% Set the axis locations
wd=0.19;
ht=[0.39 0.46];
ax_xpos=0.07:0.21:0.70;
ax_ypos=[0.58 0.07];

% Create a colour map
C=flip(colormap(pink),1);
C=[0.7*ones(64,3); repmat(C(2,:),64,1); C(2:2:256,:)];
%==========================================================================

% For each species ...
for spp=1:length(sppList)

% For each strain ...
    for s=1:length(sLab)
        disp([' ' sppList{spp} ', ' sLab{s}])
    
%==========================================================================
% COMPUTE THE SOBOL INDICES
% Create an array to store the Sobol indices
        sobol=NaN(length(parLab{spp}),length(T),nReps);

% For each replicate ...
        for r=1:nReps

% Compute the first order or total Sobol indices
            if strcmp(sensMeas,'Sobol1')==1
                sobol(:,:,r)=computeSobolIndices(T,sppListToo{spp},...
                                                 s,nSamp);
            elseif strcmp(sensMeas,'SobolT')==1
                [~,sobol(:,:,r)]=computeSobolIndices(T,sppListToo{spp},...
                                                     s,nSamp);
            end

        end

% Store the Sobol indices for each parameter
        sensInd=median(sobol,3);

% Adjust any sensitivity measures that are NaNs to so they will be grey for
% plotting
        sensInd(sensInd<0)=0;           % possible due to MC approach
        sensInd(isnan(sensInd))=-1;
%==========================================================================

%==========================================================================
% PLOT SENSITIVITY MEASURES VS TEMPERATURE
% Select the figure
        subplot('position',[ax_xpos(s) ax_ypos(spp) wd ht(spp)])
        
% Plot the sensitivity measures (note: ranges set so non-computed values
% appear as grey)
        imagesc(sensInd,[-1 1]);
        colormap(C)

% Add some dividing lines
        for j=1:length(parLab{spp})-1
            line([0.5 length(T)+0.5],(j+0.5)*[1 1],'linestyle','-',...
                 'color','k','linewidth',0.25)
        end
        
% Add a colour bar
        if s==4
            colorbar;
            ax=colorbar;
            ax.FontSize=7;
            if strcmp(sensMeas,'PRCC')==1
                ax.YLim=[-1 1];
                ax.YTick=-1:0.5:1;
                ax.YTickLabel={'-1.0','-0.5','0','0.5','1.0'};
            elseif strcmp(sensMeas(1:5),'Sobol')==1
                ax.YLim=[0 1];
                ax.YTick=0:0.2:1;
                ax.YTickLabel={'0','0.2','0.4','0.6','0.8','1.0'};
            end
            ax.Label.FontSize=7;
            if strcmp(sensMeas,'Sobol1')==1
                ax.Label.String='first-order sensitivity index';
            elseif strcmp(sensMeas,'SobolT')==1
                ax.Label.String='total sensitivity index';
            end
            os=0.5*(ht(spp)-ht(1));
            ax.Position=[ax_xpos(s)+wd+0.01 ax_ypos(spp)+os 0.02 ht(1)];
        end
        
% Make the axes look pretty
        ax=gca;
        ax.FontSize=7;
        ax.Box='on';
        ax.XTick=1:5:length(T);
        ax.XTickLabel=T(1:5:length(T));
        ax.YTick=1:length(parLab{spp});
        if s==1
            ax.YTickLabel=parLab{spp};
        else
            ax.YTickLabel={''};
        end

% Label the x-axis
        if spp==2
            xlabel('temperature ({\circ}C)')
        end

% Label the y-axis
        if s==1
            ylabel(['parameter (' sppList{spp} ')'])
        end

% Indicate the strain
        if spp==1
            title(sLab{s},'fontsize',7,'fontweight','normal')
        end
%==========================================================================

    end
end

% Print to a file
print('-dtiff','-r300',[sensMeas '_R0VsTemperatureForEHDV.tif'])

% Tidy up
close('all')
