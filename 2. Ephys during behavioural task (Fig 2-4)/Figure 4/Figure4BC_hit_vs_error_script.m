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


%% plot hit vs. error
%settings
xll = [-0.05 1.3];
yll = [-0.3 1];
xla = 'Time from stim onset (s)';
excluded_units = zeros(2,64,2,2);

clear prepData heData excluded_units_sum
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
                    prepData{task,resp}(fg,uidx,:) = nan(1,size(orgData{1,1,1},2));
                    excluded_units(task,uidx,resp,fg) = 1;
                else
                    prepData{task,resp}(fg,uidx,:) = nanmean(catData,1); %avg within cell and store
                end
            end
        end
    end
    excluded_units_sum(task,:) = squeeze(sum(sum(excluded_units(task,:,:,:),3),4));
    n(task) = sum(1-(excluded_units_sum(task,:)==4));
    
    for resp = 1:2
        subplot(2,1,resp)
        
        %save the data
        heData{task,resp} = prepData{task,resp};
        
        %plot
        for fg = 1:2 
           %average across cells and smooth the lines
            plotData_tosm(fg,:) = squeeze(nanmean(prepData{task,resp}(fg,:,:),2));
            plotData(fg,:) = smooth(plotData_tosm(fg,:),smflick);
            plotDataSem_tosm(fg,:) = squeeze(nansem_large(squeeze(prepData{task,resp}(fg,:,:)),1));
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
dpssembarfl = dpssembar';
dpsmeanbarfl = dpsmeanbar';
b = bar(1:2,dpsmeanbar,'grouped');
nBars = size(dpsmeanbarfl,1);
nGroups = size(dpsmeanbarfl,2);
groupWidth = min(0.8, nBars/(nBars + 1.5));
for i = 1:nBars
    x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);
    errorbar(x, dpsmeanbarfl(i,:),dpssembarfl(i,:), '.k');
end
xticks(1:2)
xticklabels(taskNames(1:2))
title('d-prime')
ylabel('dprime')
ylim([-1 1.4])
legend(heNames)




%% Scatter plot for figure-ground modulation hit vs. error in figure 4
%Also used this code for MP cells in supp fig

he_scatter_data = nan(2,2,2,64);
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
            he_scatter_data(task,response,fg,:) = heDataOrd(trColOut,:)';
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
        


            
            
            
           


