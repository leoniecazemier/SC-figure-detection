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

for rf = [1]
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