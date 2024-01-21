%% This script plots the amount of hits/errors/misses 

loc = 'C:\Users\leoni\Documents\NIN\Data\Optomice';
cd(loc);
logfiles = dir('*.mat');
plotit = 0;

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
yLims = {[0.45 0.95];[0.45 0.95];[0.45 0.95]};
disp_1=0;
disp_2=0;
%make colors
BlueColoravg = [0.1 0.8 1];
GreyColoravg = [0 0 0];
GreyColors = {[0.6 0.6 0.6] [0.525 0.525 0.525] [0.45 0.45 0.45] [0.375 0.375 0.375] [0.3 0.3 0.3] [0.225 0.225 0.225] [0.15 0.15 0.15] [0.075 0.075 0.075] GreyColoravg};
FGcolors = {[0 210 75]./256 [235 32 39]./256 [68 150 200]./256};
tasknames = {'Contrast','Orientation', 'Phase'};

%% 


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
%             %get p-value for performance vs. chance level
%             b.(mn).basePchance(i) = binocdf(b.(mn).baseOptoH(i), b.(mn). baseOptoH(i)+ ...
%                 b.(mn).baseOptoE(i) + b.(mn).baseOptoM(i), 0.5, 'upper') + ...
%                 binopdf(b.(mn).baseOptoH(i), b.(mn). baseOptoH(i)+ ...
%                 b.(mn).baseOptoE(i) + b.(mn).baseOptoM(i), 0.5);
        end
        %do bonferroni correction
        b.(mn).basePperf = b.(mn).basePperf .* (length(b.(mn).baseOptoH));
%         b.(mn).basePchance = b.(mn).basePchance .* (length(b.(mn).baseOptoH));
        
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
            %get p-value for performance vs. chance level
%             b.(mn).fullPchanceOri(i) = binocdf(b.(mn).fullOptoHOri(i), b.(mn). fullOptoHOri(i)+ ...
%                 b.(mn).fullOptoEOri(i) + b.(mn).fullOptoMOri(i), 0.5, 'upper') + ...
%                 binopdf(b.(mn).fullOptoHOri(i), b.(mn). fullOptoHOri(i)+ ...
%                 b.(mn).fullOptoEOri(i) + b.(mn).fullOptoMOri(i), 0.5);
        end
        %do bonferroni correction
        b.(mn).fullPperfOri = b.(mn).fullPperfOri .* (length(b.(mn).fullOptoHOri));
%         b.(mn).fullPchanceOri = b.(mn).fullPchanceOri .* (length(b.(mn).fullOptoHOri));
%         
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
            %get p-value for performance vs. chance level
%             b.(mn).fullPchanceOOP(i) = binocdf(b.(mn).fullOptoHOOP(i), b.(mn). fullOptoHOOP(i)+ ...
%                 b.(mn).fullOptoEOOP(i) + b.(mn).fullOptoMOOP(i), 0.5, 'upper') + ...
%                 binopdf(b.(mn).fullOptoHOOP(i), b.(mn). fullOptoHOOP(i)+ ...
%                 b.(mn).fullOptoEOOP(i) + b.(mn).fullOptoMOOP(i), 0.5);
        end
        %do bonferroni correction
        b.(mn).fullPperfOOP = b.(mn).fullPperfOOP .* (length(b.(mn).fullOptoHOOP));
%         b.(mn).fullPchanceOOP = b.(mn).fullPchanceOOP .* (length(b.(mn).fullOptoHOOP));
    end
    
    
end

%% summarize data

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

%% plot

fi = figure;
fi.Position = [100,100,1200,400];
colorsHEM = {[0 210 75]./256 ,[235 32 39]./256,'k'};

clear handle_mean handle_eb

for stim = 1:3
    
    for m = 1:mousetot
        mn = mice{m};
        
        if nansum(NumTrsTot(m,stim,1:12)) <trialThreshold
            NumTrsTot(m,stim,:) = NaN;
            NumHits(m,stim,:) = NaN;
            NumErrors(m,stim,:) = NaN;
            NumMisses(m,stim,:) = NaN;
            continue
        end
    end

    
    %select only latencies with at least data from 3 mice
    [r,c] = find(~isnan(squeeze(NumTrsTot(:,stim,:))));
    cols = unique(c);
    counts = histc(c, cols);
    selectLats = cols(find(counts > 2));
    usedLatsSel = myLats(selectLats);
    totals_for_lats = squeeze(NumTrsTot(:,stim,selectLats));
    hits_for_lats = squeeze(NumHits(:, stim,selectLats));
    hits_for_lats_perc = hits_for_lats ./ totals_for_lats*100;
    errors_for_lats = squeeze(NumErrors(:, stim,selectLats));
    errors_for_lats_perc = errors_for_lats ./ totals_for_lats*100;
    misses_for_lats = squeeze(NumMisses(:,stim,selectLats));
    misses_for_lats_perc = misses_for_lats ./ totals_for_lats*100;

    
    %some settings
    searchGrid.alpha = 10:1:200; %inflection point
    searchGrid.beta = -0.25:.025:0.25; %slope
    plotting_xrange = -10:200;
    plotting_yrange = [0 1];
    PF = @PAL_Logistic;
    
    %plot it
    subplot(1,3,stim)
    cla
    hold on;
    
    allreactions = cat(3,hits_for_lats,errors_for_lats,misses_for_lats);
    allreactions_perc = cat(3,hits_for_lats_perc,errors_for_lats_perc,misses_for_lats_perc);
    
    for reaction = 1:3 %hit error miss
        avg_perf = nanmean(allreactions_perc(:,:,reaction),1);
        avg_hitnum = round(avg_perf);
        avg_outofnum = ones(1,length(usedLatsSel)).*100;
        noOpPerf = avg_perf(end)/100;

        guesslim = mean(avg_hitnum(usedLatsSel<=50)./avg_outofnum(usedLatsSel<=50));
        searchGrid.lambda= (1-noOpPerf-0.05):0.005:(1-noOpPerf+0.05); %lapse rate
        searchGrid.gamma = (guesslim-0.05):0.005:(guesslim+0.05); %guessrate
        [avg_paramsValues] = PAL_PFML_BruteForceFit(usedLatsSel,avg_hitnum,avg_outofnum,searchGrid,PF);
        boon = PAL_Logistic(avg_paramsValues,plotting_xrange);

        %plot
        rPerfs = allreactions_perc(:,:,reaction);
        psems = [];
        for l = 1:size(rPerfs,2)
            psems(l) = nansem(rPerfs(:,l));
        end
        handle_mean{reaction} = plot(usedLatsSel,nanmean(allreactions_perc(:,:,reaction),1),'o', 'MarkerFaceColor',colorsHEM{reaction},'MarkerEdgeColor', 'none');
        handle_eb{reaction} = errorbar(usedLatsSel,nanmean(allreactions_perc(:,:,reaction),1), psems,'.','color',colorsHEM{reaction});
        plot(plotting_xrange,boon*100,'-', 'Color',colorsHEM{reaction})
        title(tasknames{stim})
        
    end
    
    set(gca,'xtick',[0:50:250]);
    set(gca, 'xticklabels', {'0', '50', '100', '150', '200','no'})
    ylim([0,100])
    xlabel('Laser onset (ms)');
    ylabel('Proportion of responses (%)');
    legend([handle_mean{1},handle_mean{2},handle_mean{3}],{'Hit','Error','Miss'})
    
end

