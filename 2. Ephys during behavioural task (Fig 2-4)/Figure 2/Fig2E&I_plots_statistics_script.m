%% Plots & statistics
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
    choiceRF = inputdlg(dlgQuestion2,'RF selection',1,{'[1 3]'}); %Center and edge is default
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

clear myFigs fgData hA pA stA hAreg diffValuesFromPlot
FGylims = [-0.22, 1.43; -0.22,1.43 ; -0.4, 2.2];
get_cell_ids = cell(3,3);

for rf = rfSel
    %plot preparation
    myFigs{rf} = figure(35+rf);
    set(myFigs{rf},'position', [100 500 1700 350])
    clf
    for task =1:3
        %plot preparation
        subplot(1,4,plotIdxHelp(task)+1);
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
%         switch task
%             case {1 2}
%                 ylim([-0.6 1.8]);
%             case {3}
%                 ylim([-0.9 2.7]);
%         end
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
%             selectedMiceUnique, selectedSessionsUnique,
%             selectedCellUnits, 0,nPx);   
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
    subplot(1,4,1)
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

%% Fit the onset of FGM

fittedOnset = NaN(1,3);
figure(37);
clf;
rf=1;
clear fittedLat coeff

for task=1:3
    
    tempData = fgData{rf,task};
    notCent = allRFs(whichCellSel) ~= 1;
    tempData(:,notCent,:) = [];
    
    FGM_tosm{task} = squeeze(gsubtract(tempData(1,:,:),tempData(2,:,:)));
    FGM_mean{task} = smooth(nanmean(FGM_tosm{task},1),smfvis);
    FGM_sem{task} = smooth(nansem_large(FGM_tosm{task},1),smfvis);
    
    %pixels for fitting script
    fitTimes = find(px>-0.05 & px<0.25);
    pxT = round((px(fitTimes)'+0.0005)*1000);
    [fittedLat(task),coeff(task,:)] = latencyfit_final(FGM_mean{task}(fitTimes),pxT,0);
    a = coeff(task,1); %how fast it goes to new baseline > = faster
    c = coeff(task,2); %amplitude > = higher
    d = coeff(task,3); %peak angle > = sharper
    m = coeff(task,4);  %peak time > = later
    s = coeff(task,5); %peak angle > = sharper
    fitted_FGM{task} = d * exp(m*a+0.5*s^2*a^2-a.*pxT).*normcdf(pxT,m+s^2*a,s)+c.*normcdf(pxT,m,s);
    [mf,mfix] = max(fitted_FGM{task});
    ixs = [1:mfix]; % only search before max peak
    [yL,xL]=min(abs(fitted_FGM{task}(ixs)-mf*0.33));
    
    subplot(1,3, plotIdxHelp(task))
    cla
    pxTimes = px(fitTimes);
    colFit = colorsFG{2*task -1};
    [modulsF, modulsL] = errorfill(pxTimes,FGM_mean{task}(fitTimes), FGM_sem{task}(fitTimes),colFit);
    hold on;
    plot(pxTimes,fitted_FGM{task},'k','linewidth',1.7);
    title(taskNames{task})
    axis tight
    switch task
        case {1,2}
            ylim([-0.15 0.45]);
        case 3
            ylim([-0.8 2.4]);
    end
    ylims = get(gca,'ylim');
    ylim(ylims);
    myYdiv = (ylims(2)-ylims(1))/30;
    fitAR = arrow([pxTimes(xL) fitted_FGM{task}(xL)], [pxTimes(xL) ylims(1) - myYdiv*2 ]); % Fitted inflection point
    arrow(fitAR,'EdgeColor','k','FaceColor','k');
    ylim([ylims(1)-myYdiv*2 ylims(2)]);
    text(pxTimes(xL)+0.02, ylims(1),[num2str(round(pxT(xL))) ' ms'],'color', colFit);
    zeroliney = line([min(px(fitTimes)) max(px(fitTimes))], [0 0], 'color','k','linestyle' ,'--','handleVisibility','off');
    xlabel('Time from stim onset (s)');
    ylabel('Normalized modulation');
    
    
end

%% plot hit vs. error
%settings
xll = [-0.05 1.3];
yll = [-0.3 1];
xla = 'Time from stim onset (s)';
excluded_units = zeros(2,64,2,2);

clear fgData heData excluded_units_sum
for task = 1:2
    figure(25+task);
    clf
    
    for resp = 1:2 %plot for hit and error
        clear plotData_tosm plotData plotDataSem_tosm plotDataSem
        hold on;
        
        for fg = 1:2 %figure vs. ground
            %get data
            selData = orgData(task,:,fgTypeSplit(fg,resp));
            for uidx = 1:size(selData,2)
                %concatenate desired data
                catData = selData{1,uidx};
                if size(catData,1)<mintrs
                    fgData{task,resp}(fg,uidx,:) = nan(1,size(orgData{1,1,1},2));
                    excluded_units(task,uidx,resp,fg) = 1;
                else
                    fgData{task,resp}(fg,uidx,:) = nanmean(catData,1); %avg within cell and store
                end
            end
        end
    end
    excluded_units_sum(task,:) = squeeze(sum(sum(excluded_units(task,:,:,:),3),4));
    n(task) = sum(1-(excluded_units_sum(task,:)==4));
    
    for resp = 1:2
        subplot(2,1,resp)
        
        %save the data
        heData{task,resp} = fgData{task,resp};
        
        %plot
        for fg = 1:2 
           %average across cells and smooth the lines
            plotData_tosm(fg,:) = squeeze(nanmean(fgData{task,resp}(fg,:,:),2));
            plotData(fg,:) = smooth(plotData_tosm(fg,:),smflick);
            plotDataSem_tosm(fg,:) = squeeze(nansem_large(squeeze(fgData{task,resp}(fg,:,:)),1));
            plotDataSem(fg,:) = smooth(plotDataSem_tosm(fg,:), smflick);
            
            %plot
            colidx = (task-1)*2 + fg;
            [handleFill(fg),handleLine(fg)] = errorfill(px,plotData(fg,:),plotDataSem(fg,:),colorsFG{colidx});
        end
        mlickline = line([meanRT meanRT],[-2 4] , 'color',[0.5 0.5 0.5],'linestyle' ,'--', 'linewidth',1.2, 'HandleVisibility','off');
        ylim(yll);
        if strcmpi(choiceC,'servo')
            ylim([-0.75 2.25])
        end
        xlim(xll)
        yLimits = get(gca, 'ylim');
        title([taskNames{task} ', ', heNames{resp} ',  n=' num2str(n(task))]);
        xlabel(xla);
        ylabel('Normalized firing rate');
        set(gca, 'tickdir', 'out');
        legend(handleLine,fgNames,'location','northeast')
        
    end
end



%% Check the statistics for the period between stim onset and lick
%Ori and oop together

clear heDataOrd heDataOrdCell

%get data, need to go back to orgData because I need to use lick times of
%individual trials
for task = 1:2
    for resp = 1:2 %hit/error
        for fg = 1:2
            %order as figureContra, figureIpsi, GroundContra, groundIpsi
            if fg ==1
                trIdx = resp;
            else
                trIdx = (3-resp)+2;
            end
            trIdxIn  = fgTypeSplit(fg,resp);
            trColOut = trIdx + ((task-1)*4);
            for uidx = 1:size(heData{1},2)
                if size(orgData{task,uidx,trIdxIn},1) <mintrs
                    heDataOrd(trColOut,uidx,:) = NaN;
                    heDataOrdCell{trColOut,uidx} = NaN;
                else
                    % get time between stim onset and lick
                    trGrps = unit(whichCellSel(uidx)).trialGroups;
                    pickTr = find(trGrps == (task-1)*6 + trIdxIn);
                    
                    clear tempUDat tempUDatLong
                    for trIdx = 1:length(pickTr)
                        trVal = pickTr(trIdx);
                        lickT = (unit(whichCellSel(uidx)).RTs(trVal))/1000;
                        intW = find(px>0 & px<lickT);
                        tempUDat(trIdx) = nanmean(orgData{task,uidx,trIdxIn}(trIdx,intW),2);
                        tempUDatLong(trIdx,:) = orgData{task,uidx,trIdxIn}(trIdx,:);
                    end
                    
                    heDataOrd(trColOut,uidx) = nanmean(tempUDat);
                    heDataOrdCell{trColOut,uidx} = tempUDatLong;
                end
            end
        end
    end
end

fullDataC = heDataOrd';

%make a table for mixed linear effects model
tblFR = table();
meas = 0;
taskorder = {'Orientation', 'Orientation',  'Orientation',  'Orientation', 'Phase',  'Phase',  'Phase',  'Phase'};
fgorder =   {'Figure',      'Figure',       'Ground',       'Ground',      'Figure', 'Figure', 'Ground', 'Ground'};
lickorder = {'Contra',      'Ipsi',         'Contra',       'Ipsi',        'Contra', 'Ipsi',   'Contra', 'Ipsi' };
for uidx=1:size(fullDataC,1)
    for colu = 1:size(fullDataC,2)
        if ~isnan(fullDataC(uidx,colu))
            meas = meas+1;
            Tnew = table(meas, uidx, taskorder(colu), fgorder(colu), lickorder(colu), fullDataC(uidx,colu), ...
                selectedCellMice(uidx), selectedCellSessions(uidx), ...
                'VariableNames',{'Measure','Unit','Task','FigGnd','LickSide','FiringRate','Mouse','Session'});
            tblFR = [tblFR ; Tnew];
        end
    end
end

%make models
clear lmeFR compRes
lmeFR_noint = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{1} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + Task*FigGnd + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{2} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + Task*LickSide + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{3} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + FigGnd*LickSide + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{4} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + FigGnd*LickSide + Task*FigGnd + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{5} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + Task*LickSide + Task*FigGnd + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{6} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + Task*LickSide + FigGnd*LickSide + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lmeFR{7} = fitlme(tblFR,'FiringRate ~  Task + FigGnd + LickSide + Task*FigGnd + Task*LickSide + FigGnd*LickSide + Task*LickSide*FigGnd + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');

%compare models
for i=1:length(lmeFR)
    compRes{i} = compare(lmeFR_noint, lmeFR{i},'CheckNesting', true);
    compRes_p(i) = compRes{i}.pValue(2);
end


if any(compRes_p<0.05)
    warning('Code new models for FR!')
else
    %compare the smaller model with the biggest model
    bestModel = lmeFR_noint;
end

pMainTask = bestModel.Coefficients.pValue(2);
pMainFG = bestModel.Coefficients.pValue(3);
pMainLS = bestModel.Coefficients.pValue(4);

clear bprep
bprep(1,:,:)= fullDataC(:,1:4);
bprep(2,:,:) = fullDataC(:,5:8);
for task = 1:2
  
    figure(73);
    subplot(2,1,task)
    cla
    hold on;

    bMeans = squeeze(nanmean(bprep(task,:,:),2)); %tasks together
    bMeanOrd = bMeans([1 2; 4 3]);
    bSEMs = nansem_large(squeeze(bprep(task,:,:)),1);
    bSEMOrd = bSEMs([1 2; 4 3]);
    for ll = 1:2
        mP(ll) = errorbar(bMeanOrd(ll,:),bSEMOrd(ll,:));
        mP(ll).Color = colorsFG{(task-1)*2+ll};
        mP(ll).LineWidth = 1;
    end
    xlim([0.5 2.5]);
    xticks(1:2)
    xticklabels(heNames);
    ylabel('Normalized firing rate')
    ylim([-0.1 0.4])
    title(taskNames{task})
    legend({'Figure', 'Ground'})
 
end

%% d-prime measures

%organise data as figureHit, figureError, GroundHit, groundError
dpsPrep = heDataOrdCell([1 2 4 3 5 6 8 7],:);
clear dpsInfo

for task=1:2
    for resp = 1:2 %hit error
        dpIdxOut = (task-1)*2 + resp;
        fgIdx = [resp resp+2];
        fgIdx = fgIdx + (task-1)*4;
        for uidx = 1:size(dpsPrep,2)
            %remove cells that are not included in the task
            N1 = isnan(dpsPrep{fgIdx(1),uidx}(1));
            N2 = isnan(dpsPrep{fgIdx(2),uidx}(1));
            if N1 || N2
                dpsInfo(dpIdxOut,uidx,:) = nan(1,5);
            else
                dpsInfo(dpIdxOut,uidx,:) = dprimeparts(nanmean(dpsPrep{fgIdx(1),uidx},2),... 
                    nanmean(dpsPrep{fgIdx(2),uidx},2));
            end
            %remove outliers
            if dpsInfo(dpIdxOut,uidx,1) >50 || dpsInfo(dpIdxOut,uidx,1)<-50 
                dpsInfo(dpIdxOut,uidx,:) = NaN;
            end
        end
    end
end

%prepare data for table
fullDataD = squeeze(dpsInfo(:,:,1))'; %order is Orihit, OriErr, oopHit, OOPerr
sesUnit = sesList(whichCellSel);


%make a table for mixed linear effects model
tblDp = table();
meas = 0;
taskorder = {'Orientation', 'Orientation', 'Phase', 'Phase'};
resporder = {'Hit', 'Error', 'Hit', 'Error'};
for uidx=1:size(fullDataD,1)
   sesidx = find(strcmpi(sesUnit{uidx}, sessions));
   mouseNM = sesUnit{uidx}(1:end-12);
    for colu = 1:size(fullDataD,2)
        if ~isnan(fullDataD(uidx,colu))
            meas = meas+1;
            Tnew = table(meas, uidx, sesidx, {mouseNM}, taskorder(colu), resporder(colu), fullDataD(uidx,colu), ...
                'VariableNames',{'Measure','Unit','Session', 'Mouse','Task','Response','dprime'});
            tblDp = [tblDp ; Tnew];
        end
    end
end

tblDp.Unit = categorical(tblDp.Unit);
tblDp.Session = categorical(tblDp.Session);
tblDp.Mouse = categorical(tblDp.Mouse);

fMods = {'dprime~1+Response+Task +(1|Unit)+ (1|Mouse) + (1|Session)', ...
        'dprime~1+Response*Task +(1|Unit)+ (1|Mouse) + (1|Session)'};

fmod{1} = fitlme(tblDp,fMods{1},'FitMEthod','REML','CheckHessian',1);
fmod{2} = fitlme(tblDp,fMods{2},'FitMEthod','REML','CheckHessian',1);

AICf = [fmod{1}.ModelCriterion.AIC, fmod{2}.ModelCriterion.AIC];
[minAICf, idxAICf] = min(AICf);

bestModel = fmod{idxAICf};
anModel = anova(bestModel);
bestFunc = char(bestModel.Formula);

%find contrasts for statistics of figure.
%Pull out the phase data
ixOri = find(strcmpi(tblDp.Task,'Orientation'));
ixPh = find(strcmpi(tblDp.Task,'Phase'));

%Remake table with just phase data and ensure data is categorical where
%required
TOri = table(categorical(tblDp.Unit(ixOri)),categorical(tblDp.Session(ixOri)), ... 
    tblDp.Task(ixOri),tblDp.Response(ixOri),tblDp.dprime(ixOri), tblDp.Mouse(ixOri) , ... 
    'VariableNames',{'Unit','Session','Task','Response','dprime','Mouse'});
ModelOri = fitlme(TOri,bestFunc,'FitMEthod','REML','CheckHessian',1);
POri =  ModelOri.Coefficients(2,6).pValue;

TPh = table(categorical(tblDp.Unit(ixPh)),categorical(tblDp.Session(ixPh)), ... 
    tblDp.Task(ixPh),tblDp.Response(ixPh),tblDp.dprime(ixPh), tblDp.Mouse(ixPh),... 
    'VariableNames',{'Unit','Session','Task','Response','dprime','Mouse'});
ModelPh = fitlme(TPh,bestFunc,'FitMEthod','REML','CheckHessian',1);
PPh =  ModelPh.Coefficients(2,6).pValue;


%plot the results 
figure(81);
clf
hold on;
dpsmeanbar = nanmean(fullDataD,1);
dpsmeanbar = dpsmeanbar([1 2; 3 4]);
dpssembar = nansem_large(fullDataD,1);
dpssembar = dpssembar([1 2; 3 4]);
dpssembar_flipped = dpssembar';
dpsmeanbar_flipped = dpsmeanbar';
b = bar(1:2,dpsmeanbar,'grouped');
nBars = size(dpsmeanbar_flipped,1);
nGroups = size(dpsmeanbar_flipped,2);
groupWidth = min(0.8, nBars/(nBars + 1.5));
for i = 1:nBars
    x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);
    errorbar(x, dpsmeanbar_flipped(i,:),dpssembar_flipped(i,:), '.k');
end
xticks(1:2)
xticklabels(taskNames(1:2))
title('d-prime')
ylabel('dprime')
ylim([-1 1.4])
legend(heNames)


%% scatter plots individual neurons in fig 2 & 3

%%%%%%%%%%%%%%%%%% for figure-ground modulation in figure 2 & 3

%significance windows 
contrastsigwin = {[0.07, 0.099],[],[]};
orientationsigwin = {[0.083,0.13],[],[0.08,0.103]};

fg_individual = cell(1,3); %Struct with  2*64 matrix for each task
ylims{1} = {[-1, 8],[-0.8,1.6],[-0.8,1.6]};  %rf =1
xlims{1} = {[-1, 2],[-0.8,1.6],[-0.8,1.6]};  %rf =1
ylims{3} = {[],[-0.4,2.5],[-0.4,2.5]};  %rf =3
xlims{3} = {[],[-0.4,2.5],[-0.4,2.5]};  %rf =3

for rf = [1,3]
    ifig = figure(51+rf); 
    set(ifig,'position', [100 500 1300 350])
    for task =1:3   %ori, phase, contrast  
        fg_individual{task} = nan(2,64);
        if task == 3
            sigwin = contrastsigwin{rf};
        else
            sigwin = orientationsigwin{rf};
        end
        if isempty(sigwin)
            continue
        end
        for fg = 1:2 %figure vs. ground
            %get data
            selData = orgData(task,:,fgTypeSplit(fg,:));
            for uidx = 1:size(orgData,2)
                %take out cells with wrong RF
                if allRFs(whichCellSel(uidx)) ~= rf
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
                fg_individual{task}(fg,uidx) = nanmean(nanmean(catData(:,px >= sigwin(1) & px <= sigwin(2)),2),1);
            end
        end
        subplot(1,3,plotIdxHelp(task))
        hold on; 
        colidx = (task-1)*2 + 1;
        scatter(fg_individual{task}(2,:), fg_individual{task}(1,:),20,colorsFG{colidx},'filled'); 
        line([-8, 8],[-8,8],'color','k','linestyle' ,':','handleVisibility','off')
        ylabel('Normalized Figure response')
        xlabel('Normalized Ground response')
        ylim(ylims{rf}{plotIdxHelp(task)})
        xlim(xlims{rf}{plotIdxHelp(task)})
        title(taskNames{task})
    end
end

%% Scatter plot for figure-ground modulation hit vs. error in figure 4
%Also used this code for MP cells in supp fig

he_scatter_data = nan(2,2,2,8);
ylims_he = {[-1, 1.6],[-1, 1.6]};  
xlims_he = {[-1, 1.6],[-1, 1.6]};  

plotcount = 0;

for task = 1:2  %orientation and phase
    ifig = figure(56); 
    set(ifig,'position', [100 500 1000 700])
    for response  = 1:2 %hit vs error
        plotcount = plotcount +1;
        for fg = 1:2 %figure vs ground
            if fg ==1
                trIdx = response;
            else
                trIdx = (3-response)+2;
            end
            trIdxIn  = fgTypeSplit(fg,response);
            trColOut = trIdx + ((task-1)*4);
            he_scatter_data(task,response,fg,:) = heDataOrd(trColOut,:);
        end
        
        subplot(2,2,plotcount)
        hold on;
        colidx = (task-1)*2 + 1;
        scatter(he_scatter_data(task,response,2,:), he_scatter_data(task,response,1,:),15,colorsFG{colidx});
        line([-8, 8],[-8,8],'color','k','linestyle' ,':','handleVisibility','off')
        ylabel('Normalized Figure response')
        xlabel('Normalized Ground response')
        ylim(ylims_he{task})
        xlim(xlims_he{task})
        title([taskNames{task} ', ' heNames{response}]);
    end
end
        


            
            
            
           


