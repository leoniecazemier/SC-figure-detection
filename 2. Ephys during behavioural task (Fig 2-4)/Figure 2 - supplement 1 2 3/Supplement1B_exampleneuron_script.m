%% plot example neuron

u = 54;   %neuron of choice
% load('N:\Ephys\Analysis\behTanks\recInfo');
load('C:\Users\leoni\Documents\NIN\Data\unitdatareviewed.mat')
load('C:\Users\leoni\Documents\NIN\Data\recInfo');
loc_mouse_output = 'C:\Users\leoni\Dropbox\NIN\MouseOutput\RFMap\';

trtypes = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18];  %figori grori figoop groop figgrey grgrey
colors = {[235 32 39]./256,[121 10 10]./256,[68 150 200]./256,[31 67 132]./256,[0 210 75]./256,[0 102 51]./256};
rfs = {'inside','outside', 'edge', 'unclear'};
px = -0.1995:0.001:1.7495;


%% raster plot

%get trial groups, select included trials
TG = unit(u).trialGroups;

%order trials by trialgroups
[sortedTG, TGidx] = sort(TG);
sortedTG_c = sortedTG(sortedTG>0);
TGidx_c = TGidx(sortedTG>0);

%create color categories (task, fig/gr)
colorGroups = ceil(sortedTG_c/3);

%order so that contrast task comes first
contTs = find(colorGroups >=5);
otherTs = find(colorGroups <5);
TGidx_d = [TGidx_c(contTs) TGidx_c(otherTs)];
colorGroups = [colorGroups(contTs) colorGroups(otherTs)];


%get the spikes and plot
[spR,spC] = find(unit(u).visData(TGidx_d,:));
spRCol = colorGroups(spR);
ColsP = reshape([colors{:}],[6 3]);
ColsP = ColsP([1 4 2 5 3 6],:);
ColsS = ColsP;
ColsS(:) = 0;
%also get lick times (times lick were logged from stim onset, not 200 ms
%before that)
lickTs = unit(u).RTs(TGidx_d)+200;


figure;
gscatter(spC,spR,spRCol,ColsS,'.',8)
hold on;
% scatter(lickTs,1:length(lickTs),25,[1 0.8 0],'filled');
xlabel('Time from stim onset (s)');
ylabel('Sorted Trial #');
% stimline = line([200,200],[0,max(spR)],'Color','k','linewidth',1.5,'linestyle','--');
set(gca, 'XTick', [200:500:2000], 'XTickLabel', [0:0.5:1.8]);
%select spikes from -50 ms till +250 ms
xlim([150 1450]);
set(gca, 'YDir','reverse')
set(gca ,'YTickLabel', [0:20:200])
ylim([0,max(spR)]);
grps = unique(spRCol);
for grp = grps
    thisgrp = find(spRCol == grp);
    thisR = spR(thisgrp);
    thisRU = unique(thisR);
    if ~isempty(thisRU)
%     yLimits = get(gca, 'ylim');
    xLimits = get(gca, 'xlim');
    pat(grp) = patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)  ], ... 
        [min(thisRU)-0.5 min(thisRU)-0.5 max(thisRU)+0.5  max(thisRU)+0.5], ColsP(grp,:), 'handleVisibility','off');
    alpha(pat(grp),0.25);
    set(pat(grp),'EdgeColor','none');
    end
end
box off


%% Plot Ori and phase fg modulation

clear firingrates_data
rftype = 1;
trHelp = [0 2 4];

Splits = 1; %FGM/CI
trH{1} = [1 2 3 7 8 9; 4 5 6 10 11 12];
trH{2} = [1 5 7 11; 2 4 8 10];
splName={'Figure', 'Ground'; 'Contra Lick', 'Ipsi Lick'};
colors = {[0.7 0 1],[0.2 0 0.3],[1 0.5 0],[0.5 0.15 0]};

for SPL = Splits
    figure(20+SPL)
    clf
    for ttype = 1:2
        tSelect = find(ismember(unit(u).trialGroups,trH{SPL}(ttype,:)));
        firingrates_data{ttype} = unit(u).nrDataVis(tSelect,:);
        smdatab = smooth(nanmean(firingrates_data{ttype},1),70);
        col = colors{ttype+trHelp(SPL)};
        semb = smooth(nansem_large(firingrates_data{ttype},1),70);
        [handleFill(ttype),handleLine(ttype)] = errorfill(px,smdatab,semb,col);
        %         stimline = line([0 0], [-1.5 3], 'color','k', 'linewidth',, 'HandleVisibility','off');
    end
    xticks([0:0.5:1]);
    xlabel('Time(s)');
    ylabel('Normalized firing rate');
    set(gca, 'tickdir', 'out');
    leg = splName(SPL,:);
    legend(handleLine,leg)
    ylim([-1 2]);
    xlim([-0.05 1.3])
%     title(splName{SPL});
end