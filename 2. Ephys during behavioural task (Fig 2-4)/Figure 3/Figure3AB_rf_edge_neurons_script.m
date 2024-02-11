%% plotclusteredresults3
% This script plots and analyses ephys results of Figure-ground experiments
% in Superior colliculus of mice. Spikes are spike sorted.

%% load data

% load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat')
% load('N:\Ephys\Analysis\behTanks\recInfo');
load('C:\Users\leoni\Documents\NIN\Data\unitdatareviewed.mat')
load('C:\Users\leoni\Documents\NIN\Data\recInfo');


%% settings

% This is what the trialgroups mean:
% OriFigHit    = 1;
% OriFigError  = 2;
% OriFigMiss= 3;
%
% OriGrndHit    = 4;
% OriGrndError  = 5;
% OriGrndMiss  = 6;
%
% OOPFigHit    = 7;
% OOPFigError  = 8;
% OOPFigMiss  = 9;
%
% OOPGrndHit   = 10;
% OOPGrndError  = 11;
% OOPGrndMiss = 12;
%
% GreyFigHit    = 13;
% GreyFigError  = 14;
% GreyFigMiss  = 15;
%
% GreyGrndHit    = 16;
% GreyGrndError  = 17;
% GreyGrndMiss = 18;

%rfgroups:
% inside figure = 1;
% outside figure = 2;
% edge of figure = 3;
% unclear  = 4;

% timebins
px = -0.1995:0.001:1.7495;
pxl = -1.4995:0.001:0.9995;

% Some info
nu = size(unit,2);
trtypes = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18];  %figori grori figoop groop figgrey grgrey
fgTypeSplit = [1 2 3; 4 5 6];
mintrs = 2;
sessions = unique({unit.ses});
sesList = {unit.ses};

%cell selection
oriCellSel = [unit.useOri];
oopCellSel = [unit.useOOP];
conCellSel = [unit.useOri]|[unit.useOOP];
taskCellSel = {oriCellSel oopCellSel conCellSel conCellSel};
dlgQuestion1 = 'Visual or servo motor cells?';
choiceC = questdlg(dlgQuestion1,'Cell Selection','Visual','Servo', 'Visual'); %No = default
if strcmp(choiceC,'Servo')
    cSelection = [unit.use] & [unit.MP];
    rfSel = 1:5;
else
    dlgQuestion2 = 'Which RFs to include?  1=Center;  2=Outside; 3=Edge; 4=Unclear';
    choiceRF = inputdlg(dlgQuestion2,'RF selection',1,{'[3]'}); %Center and edge is default
    rfSel = str2num(char(choiceRF{1}));
    cSelection = [unit.use] & ismember([unit.RFPOS],rfSel) & [unit.actCell];
end
cellSel = taskCellSel{4} & cSelection ;
numCellSel = sum(cellSel);
whichCellSel = find(cellSel);

%data of selected cells
allMice = {unit.mouse};
selectedCellMice = allMice(whichCellSel);
allSessions = {unit.ses};
selectedCellSessions = allSessions(whichCellSel);
selectedCellUnits = 1:64;
allRFs = [unit.RFPOS];
includedRFs = allRFs(whichCellSel);

%order and id the data of selected cells
[C,~,ib] = unique(selectedCellMice, 'stable');
selectedMiceUnique = reshape(ib, size(selectedCellMice));
[C,~,ib] = unique(selectedCellSessions, 'stable');
selectedSessionsUnique = reshape(ib, size(selectedCellSessions));

%periods of best decoding
contm = find(px>0.08 & px<=0.13);
oritm = find(px>0.08 & px<=0.13);
ooptm = find(px>0.20 & px<=0.25);
taskBestTimes = {oritm, ooptm, contm};

%for plotting
colorsFG = {[235 32 39]./256,[121 10 10]./256,[68 150 200]./256,[31 67 132]./256,[0 210 75]./256,[0 102 51]./256};
colorsP = {[190, 112, 250]/256,[95, 50, 120]/256};
rfNames = {'Inside','Outside', 'Edge', 'Unclear','All'};
fgNames = {'Figure', 'Ground'};
taskNames = {'Orientation','Phase','Contrast'};
heNames = {'Hit', 'Error'};
smfvis = 20;
smflick = 50;
plotIdxHelp = [2 3 1];

% Get mean lick time
RTS = {unit.RTs};
RTS = RTS(whichCellSel);
mRT = nan(1,length(RTS));
for cs = 1:length(RTS)
    mRT(cs) = nanmean(RTS{cs});
    minRT(cs) = min(RTS{cs});
end
meanRT = mean(mRT)/1000;


%% Prepare a variable with the desired data


%pre-allocate, task x cell x trialtype
orgData = cell(3, numCellSel,6);

for task=1:3
    for ucount = 1:numCellSel
        uReal = whichCellSel(ucount);
        
        %take out this cell for this task if performance wasn't good - don't add data
        if taskCellSel{task}(uReal) == 0
            for ttype = 1:size(orgData,3)
                orgData{task,ucount,ttype} = NaN;
            end
            continue
        end
        
        %get trial types from this task
        rowIdx = (task*2-1):task*2;
        taskTTypes = flatten(trtypes(rowIdx,:)');
        trGrps = unit(uReal).trialGroups;
        
        for tTypeIdx = 1:6
            tType = taskTTypes(tTypeIdx);
            trSel = find(trGrps == tType);  
            if ~isempty(trSel)
                orgData{task,ucount,tTypeIdx} = unit(uReal).nrDataVis(trSel,:);
            else
                orgData{task,ucount,tTypeIdx} = NaN;
            end
        end
    end
end


%% Plot FGM per task and RF

clear myFigs prepData hA pA stA hAreg diffValuesFromPlot
FGylims = [-0.22, 1.43; -0.22,1.43 ; -0.4, 2.2];
get_cell_ids = cell(3,3);

for rf = rfSel
    %plot preparation
    myFigs{rf} = figure(35+rf);
    set(myFigs{rf},'position', [100 500 1200 350])
    clf
    for task =1:2
        %plot preparation
        subplot(1,3,plotIdxHelp(task));
        cla
        hold on
        
        fgData{rf,task} = nan(2,size(orgData,2),length(px));
        
        clear plotData_tosm plotData plotDataSem_tosm plotDataSem
        for fg = 1:2 %figure vs. ground
            %get data
            selData = orgData(task,:,fgTypeSplit(fg,:));
            for uidx = 1:size(orgData,2)
                %take out cells with wrong RF
                if allRFs(whichCellSel(uidx)) ~= rf
                    fgData{rf,task}(fg,uidx,:) = NaN;
                    continue
                end
                %concatenate desired data
                catData = selData{1,uidx,1};
                if isnan(catData)
                    catData = nan(1,1950);
                end
                for ttype = 2:size(selData,3)
                    if ~isnan(selData{1,uidx,ttype})
                        catData = cat(1,catData,selData{1,uidx,ttype});
                    end
                end
                fgData{rf,task}(fg,uidx,:) = nanmean(catData,1); %avg within cell and store
            end
            
            %average across cells and smooth the lines
            plotData_tosm(fg,:) = squeeze(nanmean(fgData{rf,task}(fg,:,:),2));
            plotData(fg,:) = smooth(plotData_tosm(fg,:),smfvis);
            plotDataSem_tosm(fg,:) = squeeze(nansem_large(squeeze(fgData{rf,task}(fg,:,:)),1));
            plotDataSem(fg,:) = smooth(plotDataSem_tosm(fg,:), smfvis);
            
            %plot
            colidx = (task-1)*2 + fg;
            [handleFill(fg),handleLine(fg)] = errorfill(px,plotData(fg,:),plotDataSem(fg,:),colorsFG{colidx});
            
        end
        
        % get n
        get_cell_ids{rf,task} = whichCellSel(find(~isnan(nanmean(mean(fgData{rf,task},1),3))));
        n = length(get_cell_ids{rf,task});
        
        %log absolute differences
        diffValuesFromPlot(rf,task,:) = plotData(1,:) - plotData(2,:);
        
        %markup
        xlim([-0.05 0.25])
        title([taskNames{task} ' FGM: rf ' rfNames{rf} ', n=' num2str(n)]);
        yLimits = get(gca, 'ylim');
        legend(handleLine,fgNames);
        xlabel('Time from stim onset (s)');
        ylabel('Normalized firing rate');
        set(gca, 'tickdir', 'out');

        %find significant clusters
        ixPx = find(px>-0.05 & px<=0.25);
        nPx = px(ixPx);
        statDataFig = squeeze(fgData{rf,task}(1,:,ixPx))';
        statDataGnd = squeeze(fgData{rf,task}(2,:,ixPx))';
        
        fprintf('Running FGM statistics for rf %i; task %i. \n',rf, task)
        
        statClusters{rf,task} = ez_clusterstat_time(statDataFig,statDataGnd, 0,nPx);
%         statClusters{rf,task} = ez_clusterstat_time_LME(statDataFig,statDataGnd, ...
%             selectedMiceUnique, selectedSessionsUnique, selectedCellUnits, 0,nPx);

           % (%LME-cluster statistic is reported in paper, but
%             takes very long to run. t-test cluster is shorter, similar results)
%         
        if ~isempty(statClusters{rf,task})
            hAreg{rf,task}=regionprops(statClusters{rf,task}.map ==1,'Area','PIxelIdxList');
            clear xSigStart xSigStop
            for i=1:length(hAreg{rf,task})
                xSigStart(i) = hAreg{rf,task}(i).PixelIdxList(1);
                xSigStop(i) = hAreg{rf,task}(i).PixelIdxList(end);
                
                
                %make patch with significance
                pat = patch([nPx(xSigStart(i)) nPx(xSigStop(i)) nPx(xSigStop(i)) nPx(xSigStart(i)) ],...
                    [-0.5 -0.5 2.5 2.5], [0.5 0.5 0.5], 'handleVisibility','off');
                alpha(pat,0.2);
                set(pat,'EdgeColor','none');
            end
        end
        ylim(FGylims(3,:))
        
    end
    
    %add rf plot
    
    %preparations
    recNames = {recInfo.name};
    circx = cosd([0:1:359]);
    circy = sind([0:1:359]);
    figure(35+rf)
    subplot(1,3,1)
    hold on;
    colormap(jet)
    cmp = colormap;
    fwhm = 2*sqrt(2*log(2));
    
    for uidx = 1:size(fgData{rf,task},2)
        if allRFs(whichCellSel(uidx)) ~= rf 
            continue
        end

        uReal=whichCellSel(uidx); %unit of choice
        rfInfoCells(uidx).rfPos = rfNames{rf};
        
        %some neurons fine RFs based on visual inspection but had bad fits
        %- don't plot those fits
        if isempty(unit(uReal).params) || isempty(unit(uReal).bestFit) || (~isempty(unit(uReal).bestFit) && unit(uReal).bestFit>3)
            continue
            
        end
        
        SES = unit(uReal).ses;
        Srow = find(contains(recNames,SES));
        
        %get figure position
        figx = recInfo(Srow).figx;
        figy = recInfo(Srow).figy;
        %get RF position
        BF = unit(uReal).bestFit;
        RFx = unit(uReal).params(BF,1);
        RFy = unit(uReal).params(BF,2);
        RFwidth = unit(uReal).params(BF,3);
        
        %get relative positions
        RFrelx = RFx - figx;
        RFrely = RFy - figy;
        figXrel = figx-figx;
        figYrel = figy-figy;
        
        plot(RFrelx +circx.*fwhm.*(RFwidth/2),RFrely+circy.*fwhm.*(RFwidth/2),'color','r');
        
        xlim([-65 65])
        ylim([-60 60])
        
        %store data
        rfInfoCells(uidx).RFrelX = RFrelx;
        rfInfoCells(uidx).RFrelY = RFrely;
        rfInfoCells(uidx).RFwidth = RFwidth;
    end
    plot(figXrel+circx.*(40/2),figYrel+circy.*(40/2),'k','linewidth', 1.5)
    title('Receptive fields')
    xlabel('Rel. azimuth (vis. deg.)')
    ylabel('Rel. elevation (vis. deg.)');
end
