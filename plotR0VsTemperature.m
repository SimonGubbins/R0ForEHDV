% Matlab script to plot R0 as a function of temperature for different
% strains of epizootic haemorrhagic disease virus in cattle and deer

%==========================================================================
% PLOTTING STUFF
% List the species
sppList={'cattle'; 'deer'};
sppListToo={'Cattle'; 'Deer'};

% List the strains
sLab={'EHDV-1 (USA)',...
      'EHDV-2 (USA)',...
      'EHDV-7 (Israel)',...
      'EHDV-1 (unk.)'};
    
% Set the labels for the summary R0 measures
measLab={'max. {\it{R}}_0',...
         'temp. at max. {\it{R}}_0 ({\circ}C)',...
         'min. temp. for {{\it{R}}_0}>1 ({\circ}C)'};

% Set the plot colours
cm=colormap(lines);
col=cm([4 2 3 5],:);

% Set the width of the bars for the bar plots
bw=0.7;

% Set the y-axis limits and ticks (for the bar plots)
yAx={0:2:10, 20:26, 10:5:25};

% Set the axis locations
wd=0.21;
wd2=0.26;
ht=0.23;
ax_xpos=0.06:0.24:0.78;
ax_xpos2=0.06:0.335:0.73;
ax_ypos=[0.72 0.44 0.12];

% Specify the temperatures to plot
T=10:0.01:40;

% Specify the number of samples to draw from the joint posterior
% distribution
nReps=10;
nSamp=1000;
%==========================================================================

% For each species ...
for spp=1:length(sppList)

% For each strain ...
    for s=1:length(sLab)
        disp([' ' sppList{spp} ', ' sLab{s}])

%==========================================================================
% COMPUTE THE BASIC REPRODUCTION NUMBER
% Compute R0
        pars=[];
        for r=1:nReps
            pars=[pars; samplePosteriors(sppListToo{spp},s,nSamp)];
        end
        R0=computeR0VsTemperature(sppListToo{spp},pars,T);
%==========================================================================

%==========================================================================
% PLOT R0 VS TEMPERATURE
% Select the figure
        subplot('position',[ax_xpos(s) ax_ypos(spp) wd ht])
        
% Compute the median and 2.5 and 97.5th percentiles for R0
        pR0=prctile(R0,[50 2.5 97.5],1);

% Plot the 95% prediction interval
        patch([T flip(T,2)],[pR0(2,:) flip(pR0(3,:),2)],col(s,:),...
              'EdgeAlpha',0.5,'FaceAlpha',0.5,'linestyle','none')

% Plot the median R0
        line(T,pR0(1,:),'linestyle','-','color','k')

% Add lines indicating the threshold at R0=1
        line([min(T) max(T)],[1 1],'color','k','linestyle','--')

% Make the axes look pretty
        ax=gca;
        ax.FontSize=7;
        ax.Box='on';
        ax.XLim=[min(T) max(T)];
        ax.XTick=min(T):5:max(T);
        ax.XTickLabel=min(T):5:max(T);
        ax.YLim=[0 10];
        ax.YTick=0:2:10;
        if s==1
            ax.YTickLabel=0:2:10;
        else
            ax.YTickLabel={''};
        end

% Label the y-axis
        if s==1
            ylabel(['{\it{R}}_0 in ' sppList{spp}])
        end

% Label the x-axis
        if spp==2
            xlabel('temperature ({\circ}C)')
        end

% Indicate the strain
        if spp==1
            title(sLab{s},'fontsize',7,'fontweight','normal')
        end
%==========================================================================

%==========================================================================
% PLOT SOME SUMMARY MEASURES
% Compute the max R0, temperature at max R0 and minimum temperature for
% R0>1
        maxR0=max(R0,[],2);
        tempMaxR0=NaN(size(R0,1),1);
        minTempR0gt1=NaN(size(R0,1),1);
        for j=1:size(tempMaxR0)
            tempMaxR0(j)=mean(T(R0(j,:)==maxR0(j)));
            minT=min(T(R0(j,:)>1));
            if ~isempty(minT)
                minTempR0gt1(j)=minT;
            end
        end

% For each measure ...
        for j=1:length(measLab)

% Select the subplot
            subplot('position',[ax_xpos2(j) ax_ypos(3) wd2 ht])

% Compute the median and 2.5th and 97.5th percentiles for the measure
            if j==1
                pM=prctile(maxR0,[50 2.5 97.5],1);
            elseif j==2
                pM=prctile(tempMaxR0,[50 2.5 97.5],1);
            elseif j==3
                pM=prctile(minTempR0gt1,[50 2.5 97.5],1);
            end

% Plot a bar showing the median
            loc=length(sLab)*(spp-1)+s;
            patch(loc+0.5*bw*[-1 1 1 -1],[yAx{j}(1) yAx{j}(1) pM(1) pM(1)],...
                  col(s,:),'FaceAlpha',0.5,'EdgeColor','k','linewidth',0.5)

% Add whiskers for the 2.5th and 97.5th percentiles
            line(loc,pM(1),'linestyle','none','color','k','marker','o',...
                 'markerfacecolor','k','markersize',2)
            line(loc*[1 1],[pM(2) pM(3)],'color','k')
            
        end
%==========================================================================

    end
end

%==========================================================================
% MAKE THE SUMMARY MEASURE PLOTS LOOK PRETTY
% For each plot ...
for j=1:length(measLab)

% Select the subplot
    subplot('position',[ax_xpos2(j) ax_ypos(3) wd2 ht])

% Add a line separating cattle and deer
    line((length(sLab)+0.5)*[1 1],[0 yAx{j}(end)],'linestyle','-',...
         'color','k','linewidth',0.25)

% Add a line indicating the threshold at R0=1
    if j==1
        line([0.5 length(sppList)*length(sLab)+0.5],[1 1],'linestyle','--',...
             'color','k')
    end
    
% Set the axis properties
    ax=gca;
    ax.FontSize=7;
    ax.Box='on';
    ax.XLim=[0.5 length(sppList)*length(sLab)+0.5];
    ax.XTick=1:length(sppList)*length(sLab);
    ax.XTickLabel={''};
    ax.YLim=[yAx{j}(1) yAx{j}(end)];
    ax.YTick=yAx{j};
    ax.YTickLabel=yAx{j};

% Label the y-axis
    ylabel(measLab{j})

% Indicate the strain
    for spp=1:length(sppList)
        for s=1:length(sLab)
            yLoc=yAx{j}(1)-0.51*(yAx{j}(end)-yAx{j}(1));
            text(length(sLab)*(spp-1)+s,yLoc,sLab{s},'fontsize',5,...
                 'rotation',90,...
                 'horizontalalignment','left')
        end
    end

% Indicate the species
    yLoc=yAx{j}(1)+0.9*(yAx{j}(end)-yAx{j}(1));
    for spp=1:length(sppList)
        text(length(sLab)*(spp-1)+2.5,yLoc,sppList{spp},'fontsize',7,...
             'horizontalalignment','center')
    end

end
%==========================================================================

% Print to a file
print('-dtiff','-r300','EHDV_R0VsTemperature.tif')

% Tidy up
close('all')
clear
