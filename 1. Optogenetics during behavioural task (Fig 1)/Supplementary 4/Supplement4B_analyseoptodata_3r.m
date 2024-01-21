function [b, goodses, PS,includedMice] =  ... 
    analyseoptodata_3r(loc, bootstrapAll, nB, plotit)
%% this script takes in all log files in the given location.
% then separates them by mouse
% gets the trials where the opto tick box was on
% adds the reactions/latencies in overview variables
% different cells in reactions and optolatencies are different mice
% Adapted from a script by Lisa Kirchberger and Matt Self    

cd(loc);
logfiles = dir('*.mat');

%initialize mouse count
mouselist = {};
for f = 1:length(logfiles)
    mouselist{f} = logfiles(f).name(1:end-16);
end
mice = unique(mouselist);
mousetot = length(mice);

%some settings
b = struct();
mn = 'None';
numgoodses = 0;
goodses ={};
includedSingleMice = cell(1,3);
includedGroupMice = cell(1,3);
trialThreshold = 100;
perfThreshold= 0.5;
weightedAVG =0;
yLims = {[0.4 0.85];[0.25 0.65];[0.25 0.65]};
disp_1=0;
disp_2=0;

%make colors
BlueColoravg = [0.1 0.8 1];
GreyColoravg = [0 0 0];
GreyColors = {[0.6 0.6 0.6] [0.525 0.525 0.525] [0.45 0.45 0.45] [0.375 0.375 0.375] [0.3 0.3 0.3] [0.225 0.225 0.225] [0.15 0.15 0.15] [0.075 0.075 0.075] GreyColoravg};
FGcolors = {[0 210 75]./256 [235 32 39]./256 [68 150 200]./256};



for i=1:length(logfiles)
    ses = logfiles(i).name;
    load(ses);
    % do not use non-opto sessions
    if ~isfield(LOG, 'Optotickbox')
        continue
    elseif sum(LOG.Optotickbox) == 0
        continue
    end
    
    %check whether to include base task trials from this session
    bad = LOG.Gavepassive == 1 | LOG.Optotickbox == 0 | LOG.Drum == 1 | ...
        LOG.SemiDrum == 1 | LOG.Repeat == 1 | ~isnan(LOG.optolatency);
    
    if length(find(LOG.BGContrast ==0 & LOG.BGLuminance ==0)) >0
        gbaset = ~bad & LOG.BGLuminance == 0 & LOG.BGContrast == 0;
        hitR = sum(strcmp(LOG.Reaction, 'Hit') & strcmp(LOG.Side, 'right') & gbaset);
        hitL = sum(strcmp(LOG.Reaction, 'Hit') & strcmp(LOG.Side, 'left') & gbaset);
        errR = sum(strcmp(LOG.Reaction, 'Error') & strcmp(LOG.Side, 'right') & gbaset);
        errL = sum(strcmp(LOG.Reaction, 'Error') & strcmp(LOG.Side, 'left') & gbaset);
        perfRb = hitR/(hitR+errR);
        perfLb = hitL/(hitL+errL);
        perfGb = (hitL+hitR)/(hitR + hitL + errR + errL);
        perfDifb = abs(perfRb-perfLb);
        
        % test if both sides are significantly different from chance and not
        % different performance from each other
        pRb = binocdf(hitR, hitR+errR, 0.5, 'upper') + binopdf(hitR, hitR+errR, 0.5);
        pLb = binocdf(hitL, hitL+errL, 0.5, 'upper') + binopdf(hitL, hitL+errL, 0.5);
        
        %define whether to use the sessions
        baseThres = perfThreshold;
        useBase = perfGb >baseThres;
        if disp_1 == 0
            fprintf('Base task performance selection at %3.1f \n',baseThres);
            disp_1 = 1;
        end
        %         useBase =  perfRb > 0.6 && perfLb > 0.6 && perfGb >=0.7; %&& perfDifb <0.2;
    else
        useBase = 0;
        perfGb = [];
    end
    
    
    if length(find(LOG.BGContrast ==1 & LOG.BGLuminance ==0)) >0
        %check whether to include full task trials (ORI) from this session
        nanoop = find(isnan(LOG.OOP));
        LOG.OOP(nanoop) =0;  %OOP is not possible if mode <4 so by definition is 0 if it's logged as NAN
        gfulltOri = ~bad & LOG.BGLuminance == 0 & LOG.BGContrast == 1 & ~(isnan(LOG.OOP)|LOG.OOP);
        hitR = sum(strcmp(LOG.Reaction, 'Hit') & strcmp(LOG.Side, 'right') & gfulltOri);
        hitL = sum(strcmp(LOG.Reaction, 'Hit') & strcmp(LOG.Side, 'left') & gfulltOri);
        errR = sum(strcmp(LOG.Reaction, 'Error') & strcmp(LOG.Side, 'right') & gfulltOri);
        errL = sum(strcmp(LOG.Reaction, 'Error') & strcmp(LOG.Side, 'left') & gfulltOri);
        perfRfOri = hitR/(hitR+errR);
        perfLfOri = hitL/(hitL+errL);
        perfGfOri = (hitL+hitR)/(hitR + hitL + errR + errL);
        perfDiffOri = abs(perfRfOri-perfLfOri);
        
        % test if both sides are significantly different from chance and not
        % different performance from each other
        pRfOri = binocdf(hitR, hitR+errR, 0.5, 'upper') + binopdf(hitR, hitR+errR, 0.5);
        pLfOri = binocdf(hitL, hitL+errL, 0.5, 'upper') + binopdf(hitL, hitL+errL, 0.5);
        
        %check whether to include full task trials (OOP) from this session
        gfulltOOP = ~bad & LOG.BGLuminance == 0 & LOG.BGContrast == 1 & (LOG.OOP & ~isnan(LOG.OOP));
        hitR = sum(strcmp(LOG.Reaction, 'Hit') & strcmp(LOG.Side, 'right') & gfulltOOP);
        hitL = sum(strcmp(LOG.Reaction, 'Hit') & strcmp(LOG.Side, 'left') & gfulltOOP);
        errR = sum(strcmp(LOG.Reaction, 'Error') & strcmp(LOG.Side, 'right') & gfulltOOP);
        errL = sum(strcmp(LOG.Reaction, 'Error') & strcmp(LOG.Side, 'left') & gfulltOOP);
        perfRfOOP = hitR/(hitR+errR);
        perfLfOOP = hitL/(hitL+errL);
        perfGfOOP = (hitL+hitR)/(hitR + hitL + errR + errL);
        perfDiffOOP = abs(perfRfOOP-perfLfOOP);
        
        % test if both sides are significantly different from chance and not
        % different performance from each other
        pRfOOP = binocdf(hitR, hitR+errR, 0.5, 'upper') + binopdf(hitR, hitR+errR, 0.5);
        pLfOOP = binocdf(hitL, hitL+errL, 0.5, 'upper') + binopdf(hitL, hitL+errL, 0.5);
        
        %define whether to use the sessions
        oriThres = perfThreshold;
        useOri =  perfGfOri >=oriThres;
        if disp_2 == 0
            fprintf('Orientation task performance selection at %3.1f \n',oriThres);
        end
        %         useOri =  perfRfOri > 0.6 && perfLfOri > 0.6 && perfGfOri >=0.7; %&& perfDifb <0.2;
        oopThres = perfThreshold;
        useOOP =  perfGfOOP >=oopThres;
        if disp_2 == 0
            fprintf('Phase task performance selection at %3.1f \n',oopThres);
            disp_2 = 1;
        end
        %         useOOP =  perfRfOOP > 0.6 && perfLfOOP > 0.6 && perfGfOOP >=0.7; %&& perfDifb <0.2;
    else
        useOri = 0;
        perfGfOri = [];
        useOOP = 0;
        perfGfOOP = [];
    end
    
    if useBase || useOri || useOOP
        numgoodses = numgoodses +1;
        prepses = {ses(1:end-4); useBase; useOri; useOOP};
        goodses = [goodses [prepses]];
    end
    
    %define mouse % prepare variables
    if ~strcmp(LOG.Mouse, mn)
        mn = LOG.Mouse;
        b.(mn) = struct();
        %prepare variables
        b.(mn).optolats = [];
        b.(mn).oris = [];
        b.(mn).RTs = [];
        b.(mn).contrasts = [];
        b.(mn).bgcontrasts = [];
        b.(mn).gavepas = [];
        b.(mn).drum = [];
        b.(mn).semidrum = [];
        b.(mn).OOP = [];
        b.(mn).repeated = [];
        b.(mn).bglum = [];
        b.(mn).optoOn = [];
        b.(mn).reactions = {};
        b.(mn).sides = {};
        b.(mn).useBase = [];
        b.(mn).useOri = [];
        b.(mn).useOOP = [];
        b.(mn).propR = [];
        b.(mn).perfBase = [];
        b.(mn).perfOri = [];
        b.(mn).perfOOP = [];
    end
    
    % get behaviour data
    b.(mn).reactions = [b.(mn).reactions LOG.Reaction];
    b.(mn).optolats = [b.(mn).optolats LOG.optolatency];
    b.(mn).oris = [b.(mn).oris LOG.Orientation];
    b.(mn).sides = [b.(mn).sides LOG.Side];
    b.(mn).RTs = [b.(mn).RTs LOG.RT];
    b.(mn).contrasts = [b.(mn).contrasts LOG.FigContrast];
    b.(mn).bgcontrasts = [b.(mn).bgcontrasts LOG.BGContrast];
    b.(mn).gavepas = [b.(mn).gavepas LOG.Gavepassive];
    b.(mn).drum = [b.(mn).drum LOG.Drum];
    b.(mn).semidrum = [b.(mn).semidrum LOG.SemiDrum];
    b.(mn).OOP = [b.(mn).OOP LOG.OOP];
    b.(mn).repeated = [b.(mn).repeated LOG.Repeat];
    b.(mn).bglum = [b.(mn).bglum LOG.BGLuminance];
    b.(mn).optoOn = [b.(mn).optoOn LOG.Optotickbox];
    b.(mn).useBase = [b.(mn).useBase, repmat(useBase, [1 length(LOG.Trial)])];
    b.(mn).useOri = [b.(mn).useOri, repmat(useOri, [1 length(LOG.Trial)])];
    b.(mn).useOOP = [b.(mn).useOOP, repmat(useOOP, [1 length(LOG.Trial)])];
    b.(mn).perfBase =  [b.(mn).perfBase, repmat(perfGb, [1 length(LOG.Trial)])];
    b.(mn).perfOri =  [b.(mn).perfOri, repmat(perfGfOri, [1 length(LOG.Trial)])];
    b.(mn).perfOOP =  [b.(mn).perfOOP, repmat(perfGfOOP, [1 length(LOG.Trial)])];
    b.(mn).propR = [b.(mn).propR LOG.PropR];
    
end

%% p-values per mouse, for opto performance vs. chance level performance and non-opto performance

for m = 1:mousetot
    mn = mice{m};
    dat = b.(mn);
    
    % remove bad trials: passives, opto off,
    badtrials = dat.gavepas | ~dat.optoOn | dat.drum | dat.semidrum | dat.repeated | ~(dat.propR==0.5);
    
    %analyse base task
    basetrialsWopto = dat.contrasts == 100 & dat.bglum == 0 & dat.bgcontrasts == 0 & ~badtrials & dat.useBase;
    if sum(basetrialsWopto) >2
        [behaviourTemp] = getoptoresults_3r(dat,basetrialsWopto,mn,'Base task',plotit);
        
        b.(mn).baseOptoH = behaviourTemp{1};
        b.(mn).baseOptoE = behaviourTemp{2};
        b.(mn).baseOptoP = behaviourTemp{3};
        b.(mn).baseOptoLats = behaviourTemp{4};
        b.(mn).baseOptoM = behaviourTemp{5};
        b.(mn).basePperf = [];
        b.(mn).basePchance =[];
        for i=1:length(b.(mn).baseOptoH)-1
            %get p-value for opto performance vs non-opto performance
            if b.(mn).baseOptoP(i) < b.(mn).baseOptoP(end)
                b.(mn).basePperf(i) = binocdf(b.(mn).baseOptoH(i), b.(mn).baseOptoH(i)+ ...
                    b.(mn).baseOptoE(i) + b.(mn).baseOptoM(i), b.(mn).baseOptoP(end));
            else
                b.(mn).basePperf(i) = binocdf(b.(mn).baseOptoH(i),b.(mn).baseOptoH(i)+ ...
                    b.(mn).baseOptoE(i) + b.(mn).baseOptoM(i),b.(mn).baseOptoP(end), 'upper') +  ...
                    binopdf(b.(mn).baseOptoH(i),b.(mn).baseOptoH(i)+ ...
                    b.(mn).baseOptoE(i) + b.(mn).baseOptoM(i),b.(mn).baseOptoP(end));
            end
        end
        %do bonferroni correction
        b.(mn).basePperf = b.(mn).basePperf .* (length(b.(mn).baseOptoH));
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% full task
    
    %Ori
    fulltrialsWoptoOri = dat.contrasts == 100 & dat.bglum == 0 & dat.bgcontrasts == 1 & ...
        ~badtrials & dat.OOP == 0 & dat.useOri;
    if sum(fulltrialsWoptoOri) >2
        [behaviourTemp] = getoptoresults_3r(dat,fulltrialsWoptoOri,mn,'Full Ori',plotit);
        b.(mn).fullOptoHOri = behaviourTemp{1};
        b.(mn).fullOptoEOri = behaviourTemp{2};
        b.(mn).fullOptoPOri = behaviourTemp{3};
        b.(mn).fullOptoLatsOri = behaviourTemp{4};
        b.(mn).fullOptoMOri = behaviourTemp{5};
        b.(mn).fullPperfOri = [];
        b.(mn).fullPchanceOri =[];
        for i=1:length(b.(mn).fullOptoHOri)-1
            %get p-value for opto performance vs non-opto performance
            if b.(mn).fullOptoPOri(i) < b.(mn).fullOptoPOri(end)
                b.(mn).fullPperfOri(i) = binocdf(b.(mn).fullOptoHOri(i), b.(mn).fullOptoHOri(i)+ ...
                    b.(mn).fullOptoEOri(i) + b.(mn).fullOptoMOri(i), b.(mn).fullOptoPOri(end));
            else
                b.(mn).fullPperfOri(i) = binocdf(b.(mn).fullOptoHOri(i),b.(mn).fullOptoHOri(i)+ ...
                    b.(mn).fullOptoEOri(i) + b.(mn).fullOptoMOri(i),b.(mn).fullOptoPOri(end), 'upper') +  ...
                    binopdf(b.(mn).fullOptoHOri(i),b.(mn).fullOptoHOri(i)+ ...
                    b.(mn).fullOptoEOri(i) + b.(mn).fullOptoMOri(i),b.(mn).fullOptoPOri(end));
            end

        end
        %do bonferroni correction
        b.(mn).fullPperfOri = b.(mn).fullPperfOri .* (length(b.(mn).fullOptoHOri));
 
    end
    
    
    %OOP
    fulltrialsWoptoOOP = dat.contrasts == 100 & dat.bglum == 0 & dat.bgcontrasts == 1 ...
        & ~badtrials & dat.useOOP & dat.OOP == 1;
    if sum(fulltrialsWoptoOOP) >2
        [behaviourTemp] = getoptoresults_3r(dat,fulltrialsWoptoOOP,mn,'Full OOP',plotit);
        b.(mn).fullOptoHOOP = behaviourTemp{1};
        b.(mn).fullOptoEOOP = behaviourTemp{2};
        b.(mn).fullOptoPOOP = behaviourTemp{3};
        b.(mn).fullOptoLatsOOP = behaviourTemp{4};
        b.(mn).fullOptoMOOP = behaviourTemp{5};
        b.(mn).fullPperfOOP = [];
        b.(mn).fullPchanceOOP =[];
        for i=1:length(b.(mn).fullOptoHOOP)-1
            %get p-value for opto performance vs non-opto performance
            if b.(mn).fullOptoPOOP(i) < b.(mn).fullOptoPOOP(end)
                b.(mn).fullPperfOOP(i) = binocdf(b.(mn).fullOptoHOOP(i), b.(mn).fullOptoHOOP(i)+ ...
                    b.(mn).fullOptoEOOP(i) + b.(mn).fullOptoMOOP(i), b.(mn).fullOptoPOOP(end));
            else
                b.(mn).fullPperfOOP(i) = binocdf(b.(mn).fullOptoHOOP(i),b.(mn).fullOptoHOOP(i)+ ...
                    b.(mn).fullOptoEOOP(i) + b.(mn).fullOptoMOOP(i),b.(mn).fullOptoPOOP(end), 'upper') +  ...
                    binopdf(b.(mn).fullOptoHOOP(i),b.(mn).fullOptoHOOP(i)+ ...
                    b.(mn).fullOptoEOOP(i) + b.(mn).fullOptoMOOP(i),b.(mn).fullOptoPOOP(end));
            end
        end
        %do bonferroni correction
        b.(mn).fullPperfOOP = b.(mn).fullPperfOOP .* (length(b.(mn).fullOptoHOOP));
    end
    
    
end

%% summarize data


% initialize a matrix for storing the data
for m=1:mousetot
    bootstr_params.(mice{m}) = zeros(nB,3, 4);
    bootstr_params.(mice{m})(:) = NaN;
end

noOpPerf = nan(mousetot,3);
NumHits = nan(mousetot,3,13);
NumTrsTot = nan(mousetot,3,13);

myLats = [1000/60:1000/60:200 250];
myLatsR = [round(1000/60:1000/60:200) 250];

for stim = 1:3
    for m = 1:mousetot
        mn = mice{m};
        dat = b.(mn);
        
        if stim == 1 %base task
            if isfield(dat, 'baseOptoLats')
                useLats = find(dat.baseOptoLats <=0.25);
                StimLevels = dat.baseOptoLats(useLats) *1000;
                
                
                for LL = 1:length(useLats)
                    stimLevel = StimLevels(LL);
                    latIdx = find(myLatsR == round(stimLevel));
                    %Number of positive responses (e.g., 'yes' or 'correct' at each of the entries of 'StimLevels'
                    NumHits(m,stim,latIdx) = dat.baseOptoH(LL);
                    NumErrors(m,stim,latIdx) = dat.baseOptoE(LL);
                    NumMisses(m,stim,latIdx) = dat.baseOptoM(LL);
                    %Number of trials at each entry of 'StimLevels'
                    NumTrsTot(m,stim,latIdx) = dat.baseOptoH(LL) + dat.baseOptoE(LL) +dat.baseOptoM(LL);
                end
            else
                continue
            end
        elseif stim == 2 %Ori
            if isfield(dat, 'fullOptoLatsOri')
                useLats = find(dat.fullOptoLatsOri <=0.25);
                StimLevels = dat.fullOptoLatsOri(useLats) *1000;
                
                
                for LL = 1:length(useLats)
                    stimLevel = StimLevels(LL);
                    latIdx = find(myLatsR == round(stimLevel));
                    %Number of positive responses (e.g., 'yes' or 'correct' at each of the entries of 'StimLevels'
                    NumHits(m,stim,latIdx) = dat.fullOptoHOri(LL);
                    NumErrors(m,stim,latIdx) = dat.fullOptoEOri(LL);
                    NumMisses(m,stim,latIdx) = dat.fullOptoMOri(LL);
                    %Number of trials at each entry of 'StimLevels'
                    NumTrsTot(m,stim,latIdx) = dat.fullOptoHOri(LL) + dat.fullOptoEOri(LL) + dat.fullOptoMOri(LL);
                end
            else
                continue
            end
        elseif stim == 3 %OOP
            if isfield(dat, 'fullOptoLatsOOP')
                useLats = find(dat.fullOptoLatsOOP <=0.25);
                StimLevels = dat.fullOptoLatsOOP(useLats) *1000;
                
                
                for LL = 1:length(useLats)
                    stimLevel = StimLevels(LL);
                    latIdx = find(myLatsR == round(stimLevel));
                    %Number of positive responses (e.g., 'yes' or 'correct' at each of the entries of 'StimLevels'
                    NumHits(m,stim,latIdx) = dat.fullOptoHOOP(LL);
                    NumErrors(m,stim,latIdx) = dat.fullOptoEOOP(LL);
                    NumMisses(m,stim,latIdx) = dat.fullOptoMOOP(LL);
                    %Number of trials at each entry of 'StimLevels'
                    NumTrsTot(m,stim,latIdx) = dat.fullOptoHOOP(LL) + dat.fullOptoEOOP(LL) + dat.fullOptoMOOP(LL);
                end
            else
                continue
            end
        end
    end
end

%% do ANOVA

%test for normality
perfs = NumHits./NumTrsTot;
for task = 1:3
    perfsTask = squeeze(perfs(:,task,:));
    NumTrs = squeeze(NumTrsTot(:,task,:));
    NumTrsSel = nansum(NumTrs(:,1:12),2)>trialThreshold;
    %remove mice without data
    notnans = ~isnan(perfsTask);
    countnonans = sum(notnans,2);
    perfsTask = perfsTask(countnonans>3 &NumTrsSel,:);
    %remove latencies without data
    notnans = ~isnan(perfsTask);
    countnonans = sum(notnans,1);
    perfsTask = perfsTask(:,countnonans > 3); 

    %for anova
     nmice = size(perfsTask,1);
     nlats = size(perfsTask,2);
    
    %use anovan 1 way rm anova - get stats 
    groupsM = 1:nmice;
    groupsM = repmat(groupsM',1,nlats);
    mousenums = reshape(groupsM,[],1);
    groupsL = 1:nlats;
    groupsL = repmat(groupsL',1,nmice);
    latys = reshape(groupsL',[],1);
    [PS.avnp{task},PS.avntbl{task},PS.avnstats{task}] = anovan(reshape(perfsTask,[],1),{latys,mousenums},'random',2);
    PS.multcomp{task} = multcompare(PS.avnstats{task});
end


%% bootstrap mice together

if bootstrapAll

    origLats = [1:12] * 1000/60;
  
    ParOrNonPar = 2;
    
    Colors = {'c' 'g' 'b' 'k' 'r'};
    
    % initialize a matrix for storing the data
    bootstr_params = zeros(nB, 4);
    bootstr_params(:) = NaN;
    
    % some fitting params
    paramsFree = [1 1 1 1];  %1: free parameter, 0: fixed parameter
    PF = @PAL_Logistic;  %Alternatives: PAL_Gumbel, PAL_Weibull, PAL_CumulativeNormal, PAL_HyperbolicSecant
    searchGrid.alpha = 10:1:200; %inflection point
    searchGrid.beta = 0:.01:0.2; %slope
    
    NumTrsTotSel = NumTrsTot;
    NumHitsSel = NumHits;
    
    for stim = 1:3
        for m = 1:mousetot
             mn = mice{m};
            
            if nansum(NumTrsTot(m,stim,1:12)) <trialThreshold
                NumTrsTotSel(m,stim,:) = NaN;
                NumHitsSel(m,stim,:) = NaN;
                continue
            end
            includedGroupMice{stim,m} =  mn;
            
            minPerf(m,stim) =round(min(NumHitsSel(m,stim,:)./NumTrsTotSel(m,stim,:)),2);
            maxPerf(m,stim) =round(max(NumHitsSel(m,stim,:)./NumTrsTotSel(m,stim,:)),2);
            noOpPerf(m,stim) = NumHitsSel(m,stim,end)/NumTrsTotSel(m,stim,end);
            
            minPerf(m) =round(min(NumHitsSel(m,:)./NumTrsTotSel(m,:)),2);
            maxPerf(m) =round(max(NumHitsSel(m,:)./NumTrsTotSel(m,:)),2);
            
            
        end
        
        %don't run if there's not at least data from 2 mice for a stim
        if length(find(nansum(NumTrsTotSel(:,stim,1:12),3) > trialThreshold))<2
            continue
        end
        
        %select only latencies with at least data from 3 mice
        [r,c] = find(~isnan(squeeze(NumHitsSel(:,stim,:))));
        cols = unique(c);
        counts = histc(c, cols);
        selectLats = cols(find(counts > 2));
        selectLatsNOOP = selectLats(selectLats<13);
        if isempty(selectLats)
            continue
        end
        
        plotting_xrange = 0:201;
        plotting_yrange = [0 1];
        
        %remove latencies not used
        HitsforLats = squeeze(NumHitsSel(:, stim,selectLats));
        TotsforLats = squeeze(NumTrsTotSel(:,stim,selectLats));
        usedLatsSel = myLats(selectLats);
        %         NumPos_fittedSel = NumPos_fitted(:,selectLats);
        
        % also fit the average performance once
        if ~weightedAVG  %across mice without weighing number of trials
            avg_perf = nanmean(HitsforLats./TotsforLats,1);
            
        else %weighted avg excluding nans
            tP =  HitsforLats./TotsforLats;
            tL =  nansum(TotsforLats,2);         
            W = bsxfun(@times,~isnan(tP),tL);
            W = bsxfun(@rdivide,W,sum(W,1));
            tP(isnan(tP)) = 0;
            avg_perf = sum(tP.*W,1);
%             avg_perf = nansum((tempPerf.*temptrCount)./sum(temptrCount),1);
        end
        avg_minPerf =round(min(avg_perf),2);
        avg_maxPerf = round(max(avg_perf),2);
        avg_hitnum = round(avg_perf.*100);
        avg_outofnum = ones(1,length(usedLatsSel)).*100;
        
        %           [avg_paramsValues LL exitflag output] = PAL_PFML_Fit(usedLatsSel,avg_hitnum, ...
        %             avg_outofnum,searchGrid,paramsFree,PF, 'guessLimits', guesslimits,'lapseLimits',lapselimits);
        
        guesslim = mean(avg_hitnum(usedLatsSel<=50)./avg_outofnum(usedLatsSel<=50));
        searchGrid.lambda= (1-nanmean(noOpPerf(:, stim))-0.05):0.005:(1-nanmean(noOpPerf(:, stim))+0.05); %lapse rate
        searchGrid.gamma = (guesslim-0.05):0.005:(guesslim+0.05); %guessrate
        [avg_paramsValues ] = PAL_PFML_BruteForceFit(usedLatsSel,avg_hitnum,avg_outofnum,searchGrid,PF);
        %plot it
        if stim == 1
            figure('Position', [65   155   394   249])
        elseif stim == 2
            figure('Position', [276   162   394   249])
        elseif stim == 3
            figure('Position', [679   164   394   249])
        end
        
        %get some info for plots
        avg_p = avg_hitnum./avg_outofnum;
        myZero = find(noOpPerf == 0);
        noOpPerf(myZero) = NaN;
        mynans = ~isnan(HitsforLats);
        nanCount = sum(mynans,1);
        rPerfs = HitsforLats./TotsforLats;
        psems = [];
        for l = 1:size(rPerfs,2);
            psems(l) = nansem(rPerfs(:,l));
        end
       
        %plot it
        plot(usedLatsSel,avg_p,'o', 'MarkerFaceColor',FGcolors{stim},'MarkerEdgeColor', 'none');
        hold on
        eb = errorbar(usedLatsSel,avg_p, psems,'.','color',FGcolors{stim});
        set(gca,'xtick',[0:50:250]);
        set(gca, 'xticklabels', {'0', '50', '100', '150', '200','no'})
        boon = PAL_Logistic(avg_paramsValues,plotting_xrange);
        plot(plotting_xrange,boon,'-', 'Color',FGcolors{stim})
        xlim([0 270])
        title(['Together, stim ' num2str(stim) ', n=' num2str(max(nanCount))]);
        ylim(plotting_yrange)
        message = sprintf('Threshold estimate %f', avg_paramsValues(1));
        disp(message);
        
        %also plot individual values
        symbols = {'x', '^', 'o', 'square','diamond','*','v','pentagram', '+'};
        for mouse_plot = 1:size(rPerfs,1)
            plot(usedLatsSel,rPerfs(mouse_plot,:),symbols{mouse_plot},'color', [0.5 0.5 0.5])
        end
        
        % Bootstrapping
        % now that have all paramsValues, start the Bootstrapping
        disp('Bootstrappiiiiiiiiiiing.....');
        
        %TODO - decide how to select latencies during bootstrapping
        %    idea: keep the original latencies, and select again during
        %    each bootstrap
        
        for i = 1:nB
            
            if i == round(nB/2)
                disp('halfway done')
            end
            % generate bootstrapped data for each mouse:
            
            %prepare matrix with results
            hitnum = nan(mousetot,length(myLats));
            
            for m = 1: mousetot
                %don't bootstrap the no-opto trials
                myLatIdx = find(~isnan(NumTrsTotSel(m,stim,1:12)));
                
                for lix = 1:length(myLatIdx)
                    lat = myLatIdx(lix);
                    % make a vector with 1s for hits and 0s for errors
                    temp = [ones(NumHitsSel(m,stim,lat),1); zeros(NumTrsTotSel(m,stim,lat)-NumHitsSel(m,stim,lat),1)];
                    for k = 1 : NumTrsTotSel(m,stim,lat)
                        % generate a new distributions of 1s and 0s
                        bootstr_temp(k) = temp(randi([1 NumTrsTotSel(m,stim,lat)]));
                    end
                    %look how many hits were in the new bootstrapped data
                    hitnum(m,lat) = sum(bootstr_temp);
                    clear bootstr_temp temp
                end
                %Put in the no-opto data
                hitnum(m,13) = NumHits(m,stim,13);
            end
            
            %select only data from at least two mice - remove latencies not used
            hitnum = hitnum(:,selectLats);
            
            %take average of generated data
            if ~weightedAVG  %across mice without weighing number of trials
                avg_bootstr_perf = nanmean(hitnum./TotsforLats,1);
                
            else %weighted avg
                tP =  hitnum./TotsforLats;
                tL =  nansum(TotsforLats,2);
                W = bsxfun(@times,~isnan(tP),tL);
                W = bsxfun(@rdivide,W,sum(W,1));
                tP(isnan(tP)) = 0;
                avg_bootstr_perf = sum(tP.*W,1);
            end

            %make it into hits and outofnum just for fitting
            hitnum_all = round(avg_bootstr_perf.*100);
            outofnum_all = ones(1,length(usedLatsSel)).*100;
            
            %fit the average bootstrapped data to a new curve
            minPerf =round(min(avg_bootstr_perf),2);
            maxPerf = round(max(avg_bootstr_perf),2);
            
            %             guesslimits = [0.4 minPerf+0.05];
            %             lapselimits = [1-maxPerf-0.05 0.5];
            guesslim = mean(hitnum_all(usedLatsSel<=50)./outofnum_all(usedLatsSel<=50));
            searchGrid.gamma = (guesslim-0.05):0.005:(guesslim+0.05); %guessrate
            searchGrid.gamma = (guesslim-0.05):0.005:(guesslim+0.05); %guessrate
            
            paramsValues = PAL_PFML_BruteForceFit(usedLatsSel,hitnum_all,outofnum_all,searchGrid,PF);
            
            %save the threshold in a matrix
            bootstr_params(i,:) = paramsValues;
            
            
        end
        
        %calculate average and standard deviation:
        average_paramsSim = mean(bootstr_params,1);
        sd_paramsSim = std(bootstr_params,1);
        message = sprintf('Average Threshold Bootstrapped: %f',average_paramsSim(1));
        disp(message);
        message = sprintf('Standard Deviation of Threshold of Bootstrapping: %f',sd_paramsSim(1));
        disp(message);
        %plot it to check
        %         boon = PAL_Logistic(average_paramsSim,plotting_xrange);
        %         plot(plotting_xrange,boon,'--', 'Color', BlueColoravg)
        %         axis tight
        %         box off
        %         ylim(plotting_yrange)
        %         xlim([0 270])
        
        infPoint(stim) = average_paramsSim(1);
        infPointRounded = round(infPoint(stim));
        idxInfPoint = find(plotting_xrange == infPointRounded);
        yInfPoint = boon(idxInfPoint);
%         line([0 200] ,[0.5 0.5], 'linestyle','--','color', GreyColoravg) %chance level line
        P1spec = [infPoint(stim) yInfPoint];
        P2spec = [infPoint(stim) yLims{stim}(1)+0.025];
        AR = arrow(P1spec,P2spec); % Fitted inflection point
        arrow(AR,'EdgeColor',FGcolors{stim},'FaceColor',FGcolors{stim});  
        
        %CI
        sorted_thresholds = sort(bootstr_params(:,1));
        mean_thres(stim) = mean(sorted_thresholds);
        sem_thres(stim) = sem(sorted_thresholds);
        std_thres(stim) = std(sorted_thresholds);
        low_CI_value = round(length(sorted_thresholds)*0.025);
        high_CI_value = round(length(sorted_thresholds)*0.975);
        CI(stim,:) = [sorted_thresholds(low_CI_value) sorted_thresholds(high_CI_value)];
        message = sprintf('CI interval is %g, and %g, stdev is %2.2f', CI(stim,1), CI(stim,2), std_thres(stim));
        disp(message);
        [Hst,Pst] = lillietest(sorted_thresholds);
        normallydistributed = Pst>=0.05; 
        message2 = sprintf('Bootstrapped thresholds normally distributed: %i', normallydistributed);
        disp(message2);
        
        %         %plot CI
        %         EB = errorbar(infPoint(stim), 0.47, infPoint(stim) - CI(stim,1),CI(stim,2)- infPoint(stim),  ...
        %             'horizontal' );
        %         set(EB,'color', FGcolors{stim});
        
        % plot SEM around estimation inflection point
        EB = errorbar(infPoint(stim), yLims{stim}(1)+0.025, std_thres(stim),  'horizontal' );
        set(EB,'color', FGcolors{stim});
        
        %plot settings
        set(gca,'xtick',[0:50:250]);
        set(gca, 'xticklabels', {'0', '50', '100', '150', '200','no'})
        box off
        %         axis tight
        xlim([-10 270])
        ylim(yLims{stim})
        
    end 
    
end

%% get p-values for difference between behaviour and chance/noopto

Lats = round(100/6:100/6:200);
pChance = nan(3,5,length(Lats));
pNoop = nan(3,5,length(Lats));
mice = fieldnames(b);
mousetot = length(mice);

for m = 1:mousetot
    mn = mice{m};
    dat = b.(mn);
    %base
    if isfield(dat,'baseOptoLats')
        for ll = 1:length(dat.baseOptoLats)
            Lat = round(dat.baseOptoLats(ll)*1000);
            if Lat <= 200
                numLat = find(Lats == Lat);
%                 pChance(1,m,numLat) = dat.basePchance(ll);
                pNoop(1,m,numLat) = dat.basePperf(ll);
            end
        end
    end
    %ori
    if isfield(dat,'fullOptoLatsOri')
        for ll = 1:length(dat.fullOptoLatsOri)
            Lat = round(dat.fullOptoLatsOri(ll)*1000);
            if Lat <= 200
                numLat = find(Lats == Lat);
%                 pChance(2,m,numLat) = dat.fullPchanceOri(ll);
                pNoop(2,m,numLat) = dat.fullPperfOri(ll);
            end
        end
    end
    %oop
    if  isfield(dat,'fullOptoLatsOOP')
        for ll = 1:length(dat.fullOptoLatsOOP)
            Lat = round(dat.fullOptoLatsOOP(ll)*1000);
            if Lat <= 200
                numLat = find(Lats == Lat);
%                 pChance(3,m,numLat) = dat.fullPchanceOOP(ll);
                pNoop(3,m,numLat) = dat.fullPperfOOP(ll);
            end
        end
    end
end

for T = 1:3
    pC = squeeze(pChance(T,:,:));
    pP = squeeze(pNoop(T,:,:));
    NN = ~isnan(pC);
    NNP = ~isnan(pP);
    NS = sum(NN,1);
    NSP = sum(NNP,1);
    pickLC{T} = find(NS>2);
    pickLP{T} = find(NSP>2);
    
    pcN = pC(:,pickLC{T});
    ppN = pP(:,pickLP{T});
    
    pChanceSig{T} = pcN <0.05;
    pNoopSig{T} = ppN <0.05;
    
    pChanceSum{T} = squeeze(sum(pChanceSig{T},1));
    pNoopSum{T} = squeeze(sum(pNoopSig{T},1));
end

PS.pChanceSum = pChanceSum;
PS.pNoopSum = pNoopSum;
PS.pickLC = pickLC;
PS.pickLP = pickLP;

includedMice = {includedSingleMice includedGroupMice};
end
