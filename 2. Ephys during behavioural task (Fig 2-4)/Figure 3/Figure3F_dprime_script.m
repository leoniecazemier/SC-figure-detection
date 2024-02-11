%%
% This script is used to check out d-prime distributions from single-unit data
% It uses bootstrapping to estimate whether they are significantly
% different from zero

%% startup

load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat');

%% startup
px = -0.1995:0.001:1.7495;

%time based on decoder from ephys
csum = 0
if csum
    contm = find(px>0.00 & px<=0.16);  %peak cumulative decoding performance
    oritm = find(px>0.00 & px<=0.18); %
    ooptm = find(px>0.00 & px<=0.27); %
else
    contm = find(px>0.07 & px<=0.12);  %peak  decoding performance
    oritm = find(px>0.08 & px<=0.13); %peak decoding performance
    ooptm = find(px>0.18 & px<=0.23); %
end

myTimes = {oritm, ooptm, contm}; %which times to use for decoding
pltNames = {'Figure vs. Ground','Contra vs. Ipsi lick' };
OOPN = {'Ori','OOP','Contr'};
oriSel = [unit.useOri];
oopSel = [unit.useOOP];
conSel = [unit.useOri]|[unit.useOOP];
stSel = {oriSel oopSel conSel };

sessions = unique({unit.ses});
sesList = {unit.ses};
FFS = 1; %1 = mean 2= median
tasks = 1:3;  %1= ori 2= oop
CSS = 5 %1=rf center, 2=rf edge 3 = center, edge and outside
TBins = 1; %currently hardcoded custom timebin per task

nB = 1000;  %number of bootstraps

sessions = unique({unit.ses});
sesList = {unit.ses};

FITR = {[1 2 3];[1 5]};
GRTR = {[4 5 6];[2 4]};
decType = 1;
mintrs = 5


%% do FGM d-prime analysis

% get d-prime values

clear dprimesAllAbs dprimesAllReg dprimesBAllAbs dprimesBAllReg tracesF tracesG

for task = tasks
    selCell = stSel{task};
    addTr = 6*(task-1);
    uTime = myTimes{task};
    for cs=CSS %cell selection
        if exist('allTr','var')
            clear allTr
        end
        cellcount = 0;
        for ss = 1:length(sessions)
            ses = sessions(ss);
            %select cells - use all that qualify from this session, regardless of
            %RF location and response type.
            switch cs
                case 1 %only cells inside RF
                    usedCells = find(strcmpi(sesList,ses) & [unit.use] & selCell &  ...
                        [unit.actCell] & [unit.RFPOS] == 1);
                case 2 %only cells with RF on edge
                    usedCells = find(strcmpi(sesList,ses) & [unit.use] & selCell &  ...
                        [unit.actCell] & [unit.RFPOS] == 3);
                case 3 %only cells with RF on edge
                    usedCells = find(strcmpi(sesList,ses) & [unit.use] & selCell &  ...
                        [unit.actCell] & ismember([unit.RFPOS] ,1:3));
                case 4 %only cells with RF on edge
                    usedCells = find(strcmpi(sesList,ses) & [unit.use] & selCell &  ...
                        [unit.actCell] & ismember([unit.RFPOS] ,1:4));
                case 5 %Rf center and edge
                    usedCells = find(strcmpi(sesList,ses) & [unit.use] & selCell &  ...
                        [unit.actCell] & ismember([unit.RFPOS] ,[1 3]));
            end
            
             if isempty(usedCells) 
                continue
            end
            
            %select trials from this session
            trialVec = unit(usedCells(1)).trialGroups;
            Ftr = find(ismember(trialVec,FITR{decType}+addTr));
            Gtr = find(ismember(trialVec,GRTR{decType}+addTr));
            
            %same criterion as in decoder
            if length(Ftr)<mintrs || length(Gtr)<mintrs
                continue
            end
            
            nFtr = length(Ftr);
            nGtr = length(Gtr);
            %but only if there are enough trials in this session
            if (length(Ftr) <mintrs) || (length(Gtr) <mintrs)
                continue
            end
            %get data from specific time  trial
            
            %shuffle for bootstrap
            allTr = sort([Ftr Gtr]);
            clear FtrShuf GtrShuf
            for bb=1:nB
                FtrShuf(bb,:) = randsample(allTr,nFtr);
                GtrShuf(bb,:) = allTr(~ismember(allTr,FtrShuf(bb,:)));
            end
            
            
            for cc = 1:length(usedCells)
                %make variables for decoding
                cellcount = cellcount+1;
                uu = usedCells(cc);
                cellID{task}(cellcount) = uu;
                for mt = TBins %my time
                    uTime = myTimes{task};
                    tracesF{task,cellcount,mt} = unit(uu).nrDataVis(Ftr,uTime);
                    tracesG{task,cellcount,mt}= unit(uu).nrDataVis(Gtr,uTime);
                    
                    FF{task,cellcount,mt} = nanmean(tracesF{task,cellcount,mt},2);
                    GG{task,cellcount,mt} = nanmean(tracesG{task,cellcount,mt},2);
                    dprimesAllAbs{task,cs,mt}(cellcount) = abs(dprime(FF{task,cellcount,mt}, GG{task,cellcount,mt}));
                    dprimesAllReg{task,cs,mt}(cellcount) = dprime(FF{task,cellcount,mt}, GG{task,cellcount,mt});
                    
                    %bootstrap it!
                    for bb=1:nB
                        tracesFB = unit(uu).nrDataVis(FtrShuf(bb,:),uTime);
                        tracesGB = unit(uu).nrDataVis(GtrShuf(bb,:),uTime);
                        FFB = nanmean(tracesFB,2);
                        GGB = nanmean(tracesGB,2);
                        dprimesBAllAbs{task,cs,mt}(cellcount,bb) = abs(dprime(FFB, GGB)) ;
                        dprimesBAllReg{task,cs,mt}(cellcount,bb) = dprime(FFB, GGB) ;
                        
                    end
                end
            end
        end
    end
end

%% Plot and test the distributions of population - FGM

colors = {[235 32 39]./256,[121 10 10]./256; [68 150 200]./256,[31 67 132]./256;[0 210 75]./256,[0 102 51]./256};

mts = {'custom'}; % {'early', 'late','full vis'};
stimy = {'Orientation', 'Phase', 'Contrast'};
csels = {'RF center', 'RF edge','RF all'};
useAbsD = 0;
if useAbsD
    BootDP = dprimesBAllAbs;
    RealDP = dprimesAllAbs;
    xla = 'Absolute neuronal d';
else
    BootDP = dprimesBAllReg;
    RealDP = dprimesAllReg;
    xla = 'Neuronal d';
end

clear CI dpOrig

% bfcor = 1; %bonferroni correction
plotH = [2 3 1];
plotH_rev = [3 1 2];

for ff = FFS  %mean or median
    if ff ==1
        fu = @nanmean;
        figure(33+csum)
    else
        fu = @nanmedian;
        figure(33+csum)
    end
    for mt = 1
        for task = tasks
            tN = plotH(task);
            pln=0;
            for cs = CSS
                
%                 warning('Using absolute d-primes')
                %plot regular d-prime distribution from bootstrap
                tempDBoot = BootDP{task,cs,mt};
                tempDBoot(isinf(tempDBoot)) = NaN;
                meanPopD = fu(tempDBoot,1);

                %distribution
                minPD(task) = min(meanPopD);
                maxPD(task) = max(meanPopD);
                meanPD(task) = mean(meanPopD);
%               
                %also get distributions of dprimes from bootstrap and plot/test                
                temp = BootDP{task,cs,mt};
                temp(isinf(temp)) = NaN;
                dps{task,cs,mt} = fu(temp,1);

                %get 95% conf interval - 2 sided test
                sortedDPs= sort(dps{task,cs,mt});
                edg = [0.05 0.01 0.001]; %
                for ee = 1:length(edg)
                    edgVal = edg(ee);
                    low_CI_value = round(size(sortedDPs,2)*(0.5*edgVal));
                    high_CI_value = round(size(sortedDPs,2)*(1-(0.5*edgVal)));
                    CI{task,cs,mt,ee} = [sortedDPs(:,low_CI_value) sortedDPs(:,high_CI_value)];
                end
               
                %plot real dprime
                temp = RealDP{task,cs,mt};
                temp(isinf(temp)) = NaN;
                dpO(task) = fu(temp);
                
            end
        end
        %plot the bars
        clf
        hold on;
        CIs = squeeze(cell2mat(CI(:,CSS,1,:))); %task by low/high by significance
        CIpl=[];
        for tt=1:3
            CIpl(tt,:) = CIs(plotH_rev(tt),:,1);
        end
        bardiffs = CIpl-meanPD(plotH_rev)';
        xposs = [1:0.5:2];
        %         bar(1:3,meanPD(plotH_rev));
        myylim = [-0.5 1.75; -0.15 0.42;-0.12 0.42];
        for tt=1:3
            if tt==2
                yyaxis right
            end
            errorbar(xposs(tt), meanPD(plotH_rev(tt)),bardiffs(tt,1),bardiffs(tt,2),'.k')
            line([0.4 0.6] + tt*0.5,[dpO(plotH_rev(tt)) dpO(plotH_rev(tt))],'color',colors{plotH_rev(tt),1})
            ylim(myylim(tt,:))
            ylabel('d-prime')
        end
        line([0 2.5],[0 0],'linestyle','--','color', 'k');
        %         ylim([-0.5 1.75]);
        xlim([0.5 2.5])
        xticks(xposs)
        xticklabels(stimy(plotH_rev));
   
        
    end
end


