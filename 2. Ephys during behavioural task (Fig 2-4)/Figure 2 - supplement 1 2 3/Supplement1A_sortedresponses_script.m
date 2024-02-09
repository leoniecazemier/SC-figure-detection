%% this script makes an overview image of all used cells
% it sorts them by moment of peak activity
load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat')

%% Overview of neuron activity sorted by normalized max firing rate
% ncUse = length(unit);
whichUnits = find([unit.use]& ([unit.useOri] | [unit.useOOP]));
ncUse = length(whichUnits);

depths = [unit.depth]; 
depths = depths(whichUnits);
[sortedLocs, indSortDepth] = sort(depths);

smoothing = 1;
smF = 50;


%% get data from stim-aligned PSTH and get the order
dataPsthOrig = [];
maxlocs = [];
minlocs = [];

for i=1:length(whichUnits)
    c = whichUnits(i); 
    selection = find(ismember(unit(c).trialGroups, [1:12]));
    dat = nanmean(unit(c).convDataVis(selection,51:end-10),1);
    if smoothing
        dat = movmean(dat,smF);
    end
    dat = (dat-min(dat))./(max(dat)-min(dat));
    locmax = find(dat == max(dat));
    locmin = find(dat == min(dat));
    if length(locmax)>1
        locmax = locmax(1);
    end
    if length(locmin)>1
        locmin = locmin(1);
    end
    if isempty(locmax) %|| strcmp(unit(c).mouse, 'Cabbage') || ...
            %strcmp(unit(c).mouse, 'Flopsie')
        continue
    end
    maxlocs(i) = locmax;
    unit(c).maxloc = locmax; 
    minlocs(i) = locmin; 
    dataPsthOrig(i,:) = dat;
end

[sortedLocs, indSortLocsPsthPeak] = sort(maxlocs);
[sortedLocs, indSortLocsPsthMin] = sort(minlocs);
sortedDataPsth = dataPsthOrig(indSortLocsPsthPeak,:);
sortedDataPsthMin = dataPsthOrig(indSortLocsPsthMin,:);


%% plot the avg traces of lick and stim by each order and by depth (2x3)

%plot 
%psth by its own peak order
figure; 
imagesc(sortedDataPsth);
title('PSTHs sorted by peak')
colormap;
myc = hot; 
myc = myc(1:48,:);
colormap(gca,myc);
ylabel('Cells'); 
xlabel('Time(s)');
line([150 150] ,[0 ncUse],'color', [1 1 1]);
set(gca, 'xtick',[150:300:1850], 'xticklabel', [0:0.3:1.75], 'tickdir', 'out');