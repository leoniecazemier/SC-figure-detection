%% plot example neuron

u = 226;   %neuron of choice
%load existing unit firingrates_example
% load('N:\Ephys\Analysis\behTanks\recInfo');
% load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat');
load('C:\Users\leoni\Documents\NIN\Data\unitdatareviewed.mat')
load('C:\Users\leoni\Documents\NIN\Data\recInfo');
% loc_mouse_output = 'C:\Users\cazemier\Dropbox\NIN\MouseOutput\RFMap\';
loc_mouse_output = 'C:\Users\leoni\Dropbox\NIN\MouseOutput\RFMap\';
px = -0.1995:0.001:1.7495;


trtypes = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18];  %figori grori figoop groop figgrey grgrey
colors = {[235 32 39]./256,[121 10 10]./256,[68 150 200]./256,[31 67 132]./256,[0 210 75]./256,[0 102 51]./256};
rfs = {'inside','outside', 'edge', 'unclear'};
%% Receptive field
warning('Mind which RF (on/off/both) is the best fit!')

recNames = {recInfo.name};

SES = unit(u).ses;
Srow = find(contains(recNames,SES));

fwhm = 2*sqrt(2*log(2));

figx = recInfo(Srow).figx;
figy = recInfo(Srow).figy;
rfBl = recInfo(Srow).RFmap;
sesRF = [SES(1:end-1) num2str(rfBl)];

if u==226
    grx = figx+60;
end

load([loc_mouse_output, sesRF]);

xposs = Mat(:,2:5:end);
yposs = Mat(:,3:5:end);
xpMin = min(xposs,[],'all');
xpMax = max(xposs,[],'all');
ypMin = min(yposs,[],'all');
ypMax = max(yposs,[],'all');

gxi = linspace(xpMin,xpMax,20);
%gyi_on = linspace(min(yd_on),max(yd_on),15);
gyi = linspace(ypMin,ypMax,20);

colormap('default');
myC = colormap;
myC = [1 1 1; myC(2:end,:)];
clf;

circx = cosd([0:1:359]);
circy = sind([0:1:359]);

BF = unit(u).bestFit;  %edit here!

figure;
hold on;
rf_plot_data = squeeze(unit(u).MAPOUT(:,:,3));
FF = imagesc(gxi,gyi, rf_plot_data);
colormap(myC)

RFx = unit(u).params(BF,1);
RFy = unit(u).params(BF,2);
RFwidth = unit(u).params(BF,3);

%figure
plot(figx+circx.*(40/2),figy+circy.*(40/2),'k')

%ground
plot(grx+circx.*(40/2),figy+circy.*(40/2),'k')

%RF
plot(RFx +circx.*fwhm.*RFwidth./2,RFy+circy.*fwhm.*RFwidth./2,'r');
axis xy

axis equal
xlim([figx-22, grx+22]);
ylim([figy-30, figy+30]);

xlabel('Azimuth (vis.deg.)')
ylabel('Elevation (vis. deg.)')

ax = gca;
ax.TickDir = 'out';

%orient(gcf,'landscape');
%print('M:\Pictures firingrates_example\FGpaper\Example226_RF.pdf','-dpdf');

exRFinfo.RFxyw = [RFx RFy RFwidth];
exRFinfo.FGxy = [figx figy; grx figy];
exRFinfo.imgData = squeeze(unit(u).MAPOUT(:,:,3));
exRFinfo.imgxy = {gxi gyi};

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
ColsBG = reshape([colors{:}],[6 3]);
ColsBG = ColsBG([1 4 2 5 3 6],:);
ColsSP = ColsBG;
ColsSP(:) = 0;

figure;
gscatter(spC,spR,spRCol,ColsSP,'.',10,'.k')
xlabel('Time (s)');
ylabel('Sorted Trial #');
set(gca, 'XTick', [200:50:2000], 'XTickLabel', [0:0.05:1.5]);
%select spikes from -50 ms till +250 ms
xlim([150 450]);
set(gca, 'YDir','reverse')
set(gca ,'YTickLabel', [0:20:200])
ylim([0,max(spR)]);
%get groups
grps = unique(spRCol);
for grp = 1:length(grps)
    thisgrp = find(spRCol == grp);
    thisR = spR(thisgrp);
    thisRU = unique(thisR);
%     yLimits = get(gca, 'ylim');
    xLimits = get(gca, 'xlim');
    pat(grp) = patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)  ], ... 
        [min(thisRU)-0.5 min(thisRU)-0.5 max(thisRU)+0.5  max(thisRU)+0.5], ColsBG(grp,:), 'handleVisibility','off');
    alpha(pat(grp),0.25);
    set(pat(grp),'EdgeColor','none');
end

 
%print('M:\Pictures firingrates_example\FGpaper\Example226_Scatter.pdf','-dpdf');

%% Plot Ori and OOP fg modulation

clear firingrates_example
rftype = 1;
trHelp = [0 2 4];


for Stim = 1:3
    figure(20+Stim);
    clf
    hold on;
    for ttype = 1:2
        tSelect = find(ismember(unit(u).trialGroups, trtypes(ttype+trHelp(Stim),:)));
        firingrates_example{Stim,ttype} = unit(u).convDataVis(tSelect,:);
        smdatab = smooth(nanmean(firingrates_example{Stim,ttype},1),20);
        col = colors{ttype+trHelp(Stim)};
        semb = smooth(nansem_large(firingrates_example{Stim,ttype},1),20);
        [handleFill(ttype),handleLine(ttype)] = errorfill(px,smdatab,semb,col);
        xlim([-0.05 0.25])
        ylim([0 125]);
        yLimits = get(gca, 'ylim');
        switch Stim
            case 1
                title(['Orientation FGM example unit , RF ' rfs{unit(u).RFPOS}]);
            case 2
                title(['Phase FGM example unit , RF ' rfs{unit(u).RFPOS}]);
            case 3
                title(['Contrast FGM example unit , RF ' rfs{unit(u).RFPOS}]);
        end
        
        
    end
    
    xlabel('Time from stim onset(s)');
    ylabel('Firing rate');
    set(gca, 'tickdir', 'out');
    leg = {'Figure','Ground'};
    legend(handleLine,leg)
end

