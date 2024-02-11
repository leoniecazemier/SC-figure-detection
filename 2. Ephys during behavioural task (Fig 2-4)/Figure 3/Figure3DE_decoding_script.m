% This script is used to try and decode FG-trials from single-unit data
%In principle it will decode Figure vs. ground, but might also be used to
%decode lick side..

load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat');

%% settings
format compact
px = -0.1995:0.001:1.7495;
colors = {[235 32 39]./256,[121 10 10]./256;[68 150 200]./256,[31 67 132]./256;[0 210 75]./256,[0 102 51]./256};

sessions = unique({unit.ses});
sesList = {unit.ses};

pltNames = {'Figure vs. Ground','Contra vs. Ipsi lick' };
taskNames = {'Orientation','Phase','Contrast'};

useunit = {[unit.useOri]; [unit.useOOP]; [unit.useOri]|[unit.useOOP]};
FITR = {[1 2 3];[1 5]};
GRTR = {[4 5 6];[2 4]};
FGHE = [1, 2; 5, 4];
rfInfo = {[1],[1:3]};

taskHelp = [2 3 1];
plotHelp = [3 1 2];

myrfs = {'inside','outside', 'edge', 'unclear','all'};

slWinStarts = 0:0.010:0.25;
cellNrDataVis = {unit.nrDataVis};

verbose = 1
mintrs = 5
trNo = 10 %number of trials per trial type that are resampled
loNo = 1000 %number  of leave-one-out CVs per trialtype
rfpick = [1 3]

%% Make an artificial dataset to run the decoder on (FG Ori/OOP)
% all cells from all sessions go into one decoder. I make a large
% 'artificial' data set, that has repetitions of different (real) trials in random
% combinations. I leave one trial out and I decode this one trial based on
% the 'artifical' large dataset. We do this once for every trial to be left
% out, and then 1000x - bootstrapping to estimate variability of
% performance.

% delete(gcp('nocreate'))
% parpool(4);

clear trainMat testMat codeMat fracCor usedCells Ftr Gtr allTr meanPerf predicted p_pred p_pred_corr

for task = 1:3  %ori, oop, contrast
    for fg = 1 %2 = contra/ipsi lick
        for rfs = 1 % including all cells 
            for mt = 1:length(slWinStarts)
                
                % get trialnumbers of relevant task
                addTr = 6*(task-1);
                uTime = find(px>=slWinStarts(mt) & px< slWinStarts(mt)+0.05); %window start +50 ms window


                %for each session select correct cells and trials
                for ss = 1:length(sessions)
                    ses = sessions(ss);
                    %select cells that have RF in Figure
                    usedCells{ss} = find(strcmpi(sesList,ses) & [unit.use] & useunit{task} &  ...
                        [unit.actCell] & ismember([unit.RFPOS],rfpick));
                    
                    if isempty(usedCells{ss})
                        continue
                    end
                    %select trials from this session
                    %but only if there are enough trials in this session
                    trialVec = unit(usedCells{ss}(1)).trialGroups;
                    
                    ftr_try = find(ismember(trialVec,FITR{fg}+addTr));
                    gtr_try = find(ismember(trialVec,GRTR{fg}+addTr));

                    Ftr{ss} = ftr_try;
                    Gtr{ss} = gtr_try;
                    
                    if (length(ftr_try) <mintrs) || (length(gtr_try) <mintrs)
                        usedCells{ss}=[];
                        continue
                    end
                    
                    
                end
                
                %prep trial selection for CV
                tname = taskNames(task);
                Ftr_t = Ftr(:);
                Gtr_t = Gtr(:);
                predictTemp =[];
                
                for leaveOut = 1:loNo*2
                    
                    % monitor the script
                    if leaveOut ==1
                        btic = tic;
                    end
                    if verbose
                        %                             fprintf('LeaveOut %s %i, \n', char(tname), leaveOut);
                    end
                    if leaveOut==1
                        fprintf('Started decoding for mintrs %i, task %i, fg %i, time %i, RFs %i.  \n',mintrs, task, fg, mt, rfs);
                    end
                  
                    %create training/testing sets
                    cellcount = 0;
                    trainMat = [];
                    testMat = [];
                    if leaveOut == 1
                        unit_used{task} = [];
                    end
                    for ss=1:length(Ftr_t)
                        if isempty(usedCells{ss})
                            continue
                        end
                        
                        %leave out one trial (balanced)
                        for cc = 1:length(usedCells{ss})
                            uu = usedCells{ss}(cc);
                            if leaveOut == 1
                                unit_used{task} = [unit_used{task} uu];
                            end
                            
                            allTrPlain = [Ftr_t{ss} Gtr_t{ss}];
                            
                            if leaveOut<loNo+1
                                leftOut = randsample(Ftr_t{ss},1,0);
                                Ftr_temp = Ftr_t{ss}(Ftr_t{ss} ~= leftOut);
                                Gtr_temp = Gtr_t{ss};
                            else
                                leftOut = randsample(Gtr_t{ss},1,0);
                                Gtr_temp = Gtr_t{ss}(Gtr_t{ss} ~= leftOut);
                                Ftr_temp = Ftr_t{ss};
                            end
                            
                            % create surrogate dataset without CV trial
                            allTr_temp = [Ftr_temp Gtr_temp];
                            
                            %select trials randomly
                            FtrUsed = randsample(Ftr_temp, trNo,1);
                            GtrUsed = randsample(Gtr_temp, trNo,1);
                            %
                            
                            %Get data from specific time  trial
                            if cc==1
                                trainMatTemp = [];
                                testMatTemp = [];
                                if leaveOut==1 && leaveOut==1
                                    mattic = tic;
                                end
                            end
                            
                            %make the training and test data
                            allTr = [FtrUsed GtrUsed];
                            useTraces = cellNrDataVis{uu}(allTr,uTime);
                            
                            useVals = nanmean(useTraces,2);
                            trainMatTemp(cc,:) = useVals;
                            
                            testValFull = cellNrDataVis{uu}(leftOut,uTime);
                            testVal = nanmean(testValFull);
                            testMatTemp(cc) = testVal;
                            
                        end
                        trainMat = [trainMat; trainMatTemp];
                        testMat = [testMat testMatTemp];
                    end
                    
                    if isempty(trainMat)
                        predictTemp = [predictTemp NaN];
                        continue
                    end
                    
                    %monitor script
                    if leaveOut==1 && verbose
                        mattoc = toc(mattic);
                        fprintf('Building one trainMat takes %5.1f seconds. \n', mattoc);
                    end
                    
                    %give answers to training data
                    codeMatTrs = [repmat(1,trNo,1); repmat(2, trNo, 1)];
                    codeMatLo = [repmat(1,loNo,1); repmat(2, loNo, 1)];
                    
                    %decode
                    %prepare variables
                    trainingData =  [trainMat ; codeMatTrs'];
                    
                    %build model
                    if leaveOut ==1 && verbose
                        modeltic = tic;
                        [trainedClassifier, validationAccuracy] = trainClassifier_linSVM(trainingData,5);
                        modeltoc = toc(modeltic);
                        fprintf('Training the model takes %5.1f seconds. \n', modeltoc);
                    else
                        [trainedClassifier, validationAccuracy] = trainClassifier_linSVM(trainingData,5);
                    end
                    
                    
                    %test the model on the trial that was left out
                    if leaveOut==1 && verbose
                        predtic = tic;
                        predictTemp = [predictTemp trainedClassifier.predictFcn(testMat')];
                        predtoc = toc(predtic);
                        
                        fprintf('Predicting the test trial takes %5.1f seconds.\n', predtoc);
                    else
                        predictTemp = [predictTemp trainedClassifier.predictFcn(testMat')];
                    end
                    if leaveOut==1 && verbose
                        btoc = toc(btic);
                        fprintf('One leaveOut takes %5.1f seconds. \n', btoc);
                    end
                    
                    %save the weights
                    ModelW{task,mt,leaveOut} = trainedClassifier.ClassificationSVM.Beta;
                end
                predicted{rfs,fg,mt,task} = predictTemp;
                %calculate performance; use summation in log-domain to avoid numerical errors
                try
                    notnans = find(~isnan(predictTemp));
                    lno = length(notnans);
                    fracCor = sum(predictTemp(notnans)'== codeMatLo(notnans))/(length(notnans));
                    numCor = length(find(predictTemp(notnans)'== codeMatLo (notnans)));
                    if numCor<(lno/2)
                        p_pred(rfs,fg,mt,task) = binocdf(numCor, lno, 0.5);
                    else
                        p_pred(rfs,fg,mt,task) = binocdf(numCor, lno, 0.5, 'upper') +  binopdf(numCor, lno, 0.5);
                    end
                    meanPerf(rfs,fg,mt, task) = fracCor;
                catch
                    fprintf('could not get fraccor for mintrs %i, task %i, fg %i, time %i and RFs %i.  \n',mintrs, task, fg, mt, rfs);
                    p_pred(rfs,fg,mt,task) = NaN;
                    meanPerf(rfs,fg,mt,task) = NaN;
                end
                
                
            end
            %multiple comparison correction bonferroni-holm
            predTemp = squeeze(p_pred(rfs,fg,:,task));
            predNew = bonf_holm(predTemp);
            p_pred_corr(rfs,fg,:,task) = predNew;
        end
    end
    
end



%save the workspace (or the important part)
save(['C:\Users\leoni\Documents\NIN\Data\Decoding\decoderoutput' date '.mat'], ...
    'predicted', 'p_pred', 'p_pred_corr', 'meanPerf','unit_used','ModelW')


%% plot the interesting data

xps = 0.025:0.01:0.275;

%get maximum decoding performance bin for each task
clear rfMax tbMax tbMaxStart maxPerfVal
for fg=1
    for task = 1:3
        %get maximum performance bin
        maxPerfVal(fg,task) = max(meanPerf(1,fg,:,task),[],'all');
        meanPerfSmall = squeeze(meanPerf(1,fg,:,task));
        [tbMax(fg, task)] = find(meanPerfSmall == maxPerfVal(fg,task),1);
        tbMaxStart(fg, task) = slWinStarts(tbMax(fg,task));
    end
end

centTime = xps(tbMax(task));


for fg = 1
    figure('name', pltNames{fg})
    for task =1:3
        subplot(1,3,taskHelp(task));
        cla;
        hold on;
        %get data
        selData = squeeze(meanPerf(:,fg,:,task));
        plot(xps, selData*100, 'color', colors{task,1}, 'marker','o', ...
            'markerfacecolor',colors{task,1}, 'linewidth',1.5,'markersize', 4)
        %         plot(xps, selData(2,:), 'color', colors{task,2},'marker','^',...
        %             'markerfacecolor',colors{task,2}, 'linewidth',1.5,'markersize', 4)
        
        ylim([30 100]);
        yLimits = get(gca, 'ylim');
        
        %also make significance patches
        hAd.(taskNames{task})=regionprops(p_pred_corr(1,fg,:,task) <0.05,'PIxelIdxList','Area');
        hSizes = [hAd.(taskNames{task}).Area];
        %         hAreg(hSizes<10) = [];
        
        if ~isempty(hAd.(taskNames{task}))
            clear xSigStart xSigStop
            for i=1:length(hAd.(taskNames{task}))
                xSigStart(i) = hAd.(taskNames{task})(i).PixelIdxList(1);
                xSigStop(i) = hAd.(taskNames{task})(i).PixelIdxList(end);
                
                
                %make patch with significance
                pat = patch([xps(xSigStart(i))-0.005 xps(xSigStop(i))+0.005 xps(xSigStop(i))+0.005 xps(xSigStart(i))-0.005 ],...
                    [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.5 0.5 0.5], 'handleVisibility','off');
                alpha(pat,0.2);
                set(pat,'EdgeColor','none');
            end
        end
        
        xlim([0 0.25])
        xticks(0:0.05:0.25);
        %         xtickangle(35)
        %         legend({'RF center', 'RF all'},'location','south');
        xlimits = get(gca, 'xlim');
        line(xlimits, [50 50], 'color','k','linestyle','--', 'HandleVisibility','off')
        
        title(taskNames{task});
        ylabel('Decoding performance (%)')
        xlabel('Time from stim onset (s)')
    end
end

%% Retrieve weights and plot 
RFPOS = [unit.RFPOS];
%get sessions
clear Weights Wmean Wsem Wall

for task= 1:2
    %summarize data
    bestT = tbMaxStart(1,task);
    bestTidx = find(slWinStarts ==bestT);
    units = unit_used{task};
    RFtask = RFPOS(units);
    for LO = 1:size(ModelW,3)
        Weights{task}(LO,:) = ModelW{task,bestTidx,LO};
    end
    %get mean
    Weights_mean = -mean(Weights{task},1);
    for rf = 1:length(rfpick)
        Wmean(task,rf) = mean(Weights_mean(RFtask == rfpick(rf)));
        Wsem(task,rf) = sem(Weights_mean(RFtask == rfpick(rf)));
        Wall{task}{rfpick(rf)} = Weights_mean(RFtask == rfpick(rf));
    end
    
end


%plot it
figure(46)
clf
hold on;
WmeanOrd = Wmean';
WsemOrd = Wsem';
BAR = bar(1:2, WmeanOrd', 'grouped');
nBars = size(WmeanOrd,1);
nGroups = size(WmeanOrd,2);
groupWidth = min(0.8, nBars/(nBars + 1.5));
for i = 1:nBars
    x = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);
    errorbar(x, WmeanOrd(i,:), WsemOrd(i,:), '.k');
end
xticks(1:2)
xticklabels(taskNames)
title('Model weights')
ylabel('Relative weight')
ylim([-0.1 0.5])
legend(rfNames{[1 3]})


%% weights statistics
%cell selection
oriCellSel = [unit.useOri];
oopCellSel = [unit.useOOP];
conCellSel = [unit.useOri]|[unit.useOOP];
taskCellSel = {oriCellSel oopCellSel conCellSel conCellSel};
rfSel = [1 3];
cSelection = [unit.use] & ismember([unit.RFPOS],rfSel) & [unit.actCell];
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

%
excluded_units_lme = {[34 35 44],[21 22 23 24 25 26 27 39 40,42 43 45 46 47 48 49 50 51]};
included_units_sum = repmat([1],[2 64]);
included_units_sum(1,excluded_units_lme{1}) = 0;
included_units_sum(2,excluded_units_lme{2}) = 0;
weights_mean = nan(2,64);

for task = 1:2
    included_units_lme = find(included_units_sum(task,:));
    for uidx = 1:length(included_units_lme)
        unit_real = included_units_lme(uidx);
        weights_mean(task,unit_real) = nanmean(Weights{task}(:,uidx));
    end
end

%prepare some data
RF_names = {'Center',[],'Edge'};
includedRFs_lme = RF_names(includedRFs);
tasknames = {'Orientation','Phase'};

%make a table for mixed linear effects model
tbl = table();
meas = 0;
for uidx=1:size(weights_mean,2)
    for task = 1:2
        if ~isnan(weights_mean(task,uidx))
            meas = meas+1;
            Tnew = table(meas, uidx, tasknames(task), includedRFs_lme(uidx), weights_mean(task,uidx), ... 
                selectedCellMice(uidx), selectedCellSessions(uidx), ...                
                'VariableNames',{'Measure','Unit','Task','RF','Modelweight','Mouse','Session'});
            tbl = [tbl ; Tnew];
        end
    end
end

lme{1} = fitlme(tbl,'Modelweight ~ Task + RF + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');
lme{2} = fitlme(tbl,'Modelweight ~ Task + RF + Task * RF + (1|Unit) + (1|Session) + (1|Mouse)','FitMEthod','REML');

for i=1:length(lme)
    BIC(i) = lme{i}.ModelCriterion.BIC;
end

%compare the smaller model with the biggest model
[~, bestModel_i] = min(BIC);
bestModel = lme{bestModel_i};

% pMainTask = bestModel.Coefficients.pValue(2);
% pMainRF = bestModel.Coefficients.pValue(3);

% post-hoc
ixOri = find(strcmpi(tbl.Task,'Orientation'));
ixPh = find(strcmpi(tbl.Task,'Phase'));

%Remake table with just phase data and ensure data is categorical where
%required
TOri = table(tbl.Measure(ixOri),tbl.Unit(ixOri),tbl.Task(ixOri),tbl.RF(ixOri),tbl.Modelweight(ixOri), ... 
    'VariableNames',{'Measure','Unit','Task','RF','Modelweight'});
ModelOri = fitlme(TOri,'Modelweight ~  Task + RF + (1|Unit)','FitMEthod','REML');
% POri =  ModelOri.Coefficients(2,6).pValue;

TPh = table(tbl.Measure(ixPh),tbl.Unit(ixPh),tbl.Task(ixPh),tbl.RF(ixPh),tbl.Modelweight(ixPh), ... 
    'VariableNames',{'Measure','Unit','Task','RF','Modelweight'});
ModelPh = fitlme(TPh,'Modelweight ~  Task + RF + (1|Unit)','FitMEthod','REML');
% POri =  ModelOri.Coefficients(2,6).pValue;













