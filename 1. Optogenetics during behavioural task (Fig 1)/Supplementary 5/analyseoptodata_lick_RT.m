function [b, goodses, PS] = ...
    analyseoptodata_lick_RT(loc, plotit);
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
trialThreshold = 100;
perfThreshold= 0.5;
disp_1=0;
disp_2=0;
%make colors
GreyColoravg = [0 0 0];
colorsFG = {[0 210 75]./256,[0 102 51]./256, [235 32 39]./256,[121 10 10]./256,[68 150 200]./256,[31 67 132]./256};

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
        b.(mn).RTrightVec =[];
        b.(mn).RTleftVec = [];
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
    b.(mn).RTrightVec = [b.(mn).RTrightVec LOG.RTrightVec];
    b.(mn).RTleftVec = [b.(mn).RTleftVec LOG.RTleftVec];
   
end

%% get RTs per mouse

for m = 1:mousetot
    mn = mice{m};
    dat = b.(mn);
    
    % remove bad trials: passives, opto off,
    badtrials = dat.gavepas | ~dat.optoOn | dat.drum | dat.semidrum | dat.repeated | ~(dat.propR==0.5);
    
    %analyse base task
    basetrialsWopto = find(dat.contrasts == 100 & dat.bglum == 0 & dat.bgcontrasts == 0 & ~badtrials & dat.useBase);
    if sum(basetrialsWopto) >2
        b.(mn).baseRTs = dat.RTs(basetrialsWopto);
        b.(mn).baseReactions = dat.reactions(basetrialsWopto);
        b.(mn).baseSides = dat.sides(basetrialsWopto);
        b.(mn).baseLats = dat.optolats(basetrialsWopto) * 1000;
        b.(mn).baseLats(isnan(b.(mn).baseLats)) = 250;
        b.(mn).baseRTrightVec = dat.RTrightVec(basetrialsWopto);
        b.(mn).baseRTleftVec  = dat.RTleftVec(basetrialsWopto);
       
    end
    
    
    %%%%%%%% full task
    
    %Ori
    fulltrialsWoptoOri = find(dat.contrasts == 100 & dat.bglum == 0 & dat.bgcontrasts == 1 & ...
        ~badtrials & dat.OOP == 0 & dat.useOri);
    if sum(fulltrialsWoptoOri) >2
        b.(mn).oriRTs = dat.RTs(fulltrialsWoptoOri);
        b.(mn).oriReactions = dat.reactions(fulltrialsWoptoOri);
        b.(mn).oriSides = dat.sides(fulltrialsWoptoOri);
        b.(mn).oriLats = dat.optolats(fulltrialsWoptoOri) *1000;
        b.(mn).oriLats(isnan(b.(mn).oriLats)) = 250;
        b.(mn).oriRTrightVec = dat.RTrightVec(fulltrialsWoptoOri);
        b.(mn).oriRTleftVec  = dat.RTleftVec(fulltrialsWoptoOri);
    end
    
    
    %OOP
    fulltrialsWoptoOOP = find(dat.contrasts == 100 & dat.bglum == 0 & dat.bgcontrasts == 1 ...
        & ~badtrials & dat.useOOP & dat.OOP == 1);
    if sum(fulltrialsWoptoOOP) >2
        b.(mn).oopRTs = dat.RTs(fulltrialsWoptoOOP);
        b.(mn).oopReactions = dat.reactions(fulltrialsWoptoOOP);
        b.(mn).oopSides = dat.sides(fulltrialsWoptoOOP);
        b.(mn).oopLats = dat.optolats(fulltrialsWoptoOOP)*1000;
        b.(mn).oopLats(isnan(b.(mn).oopLats)) = 250;
        b.(mn).oopRTrightVec = dat.RTrightVec(fulltrialsWoptoOOP);
        b.(mn).oopRTleftVec  = dat.RTleftVec(fulltrialsWoptoOOP);
    end
    
    
end

%% summarize data


% initialize a matrix for storing the data
myLats = [1000/60:1000/60:200 250];
myLatsR = [round(1000/60:1000/60:200) 250];
reactions = {'Hit', 'Error'};
stims = {'Contrast', 'Orientation', 'Phase'};
clear meanRTsHE meanRTsAll

for stim = 1:3
    for m = 1:mousetot
        mn = mice{m};
        dat = b.(mn);
        
        if stim == 1 %base task
            if isfield(dat, 'baseLats')
                if sum((strcmp(dat.baseReactions,'Hit') | strcmp(dat.baseReactions,'Error')) & dat.baseLats <250) >100
                    for LL = 1:length(myLatsR)
                        stimLevel = myLatsR(LL);
                        trials = find(round(dat.baseLats) == round(stimLevel) & (strcmp(dat.baseReactions, reactions(1)) | strcmp(dat.baseReactions, reactions(2))));
                        meanRTsAll{stim}(m,LL) = mean(dat.baseRTs(trials));
                        
                        for reaction = 1:length(reactions)
                            trials = find(round(dat.baseLats) == round(stimLevel) & strcmp(dat.baseReactions, reactions(reaction)));
                            meanRTsHE{stim}(m,LL,reaction) = mean(dat.baseRTs(trials));
                        end
                    end
                else
                    meanRTsHE{stim}(m,:,:) = NaN;
                    meanRTsAll{stim}(m,:) = NaN;
                end
            else
                meanRTsHE{stim}(m,:,:) = NaN;
                meanRTsAll{stim}(m,:) = NaN;
            end
        elseif stim == 2 %Ori
            if isfield(dat, 'oriLats')
                if sum((strcmp(dat.oriReactions,'Hit') | strcmp(dat.oriReactions,'Error')) & dat.oriLats <250) >100
                    for LL = 1:length(myLatsR)
                        stimLevel = myLatsR(LL);
                        trials = find(round(dat.oriLats) == round(stimLevel) & (strcmp(dat.oriReactions, reactions(1)) | strcmp(dat.oriReactions, reactions(2))));
                        meanRTsAll{stim}(m,LL) = mean(dat.oriRTs(trials));

                        for reaction = 1:length(reactions)
                            trials = find(round(dat.oriLats) == round(stimLevel) & strcmp(dat.oriReactions, reactions(reaction)));
                            meanRTsHE{stim}(m,LL,reaction) = mean(dat.oriRTs(trials));
                        end
                    end
                else
                    meanRTsHE{stim}(m,:,:) = NaN;
                    meanRTsAll{stim}(m,:) = NaN;
                end
            else
                meanRTsHE{stim}(m,:,:) = NaN;
                meanRTsAll{stim}(m,:) = NaN;
            end
        elseif stim == 3 %OOP
            if isfield(dat, 'oopLats')
                if sum((strcmp(dat.oopReactions,'Hit') | strcmp(dat.oopReactions,'Error')) & dat.oopLats <250) >100
                    for LL = 1:length(myLatsR)
                        stimLevel = myLatsR(LL);
                        trials = find(round(dat.oopLats) == round(stimLevel) & (strcmp(dat.oopReactions, reactions(1)) | strcmp(dat.oopReactions, reactions(2))));
                        meanRTsAll{stim}(m,LL) = mean(dat.oopRTs(trials));

                        for reaction = 1:length(reactions)
                            trials = find(round(dat.oopLats) == round(stimLevel) & strcmp(dat.oopReactions, reactions(reaction)));
                            meanRTsHE{stim}(m,LL,reaction) = mean(dat.oopRTs(trials));
                        end
                    end
                else
                    meanRTsHE{stim}(m,:,:) = NaN;
                    meanRTsAll{stim}(m,:) = NaN;
                end
            else
                meanRTsHE{stim}(m,:,:) = NaN;
                meanRTsAll{stim}(m,:) = NaN;
            end
        end
    end
    meanRTsHE{stim}(meanRTsHE{stim} == 0) = NaN;
    meanRTsAll{stim}(meanRTsAll{stim} == 0) = NaN;
end

%% Plot RTs
for stim=1:3
    
    %Plot hit and error together
    figure(stim+103)
    clf
    color = 2*(stim-1)+1;
    stimRTs = meanRTsAll{stim};
    pickLats = find(sum(~isnan(squeeze(stimRTs(:,:,1))),1)>3);
    hold on;
    ph = plot(myLatsR(pickLats),nanmean(stimRTs(:,pickLats),1),'o', 'MarkerFaceColor',colorsFG{color},'MarkerEdgeColor', 'none');
    ebh = errorbar(myLatsR(pickLats),nanmean(stimRTs(:,pickLats),1), nansem_large(stimRTs(:,pickLats),1),'.','color',colorsFG{color});
    set(gca,'xtick',[0:50:250]);
    set(gca, 'xticklabels', {'0', '50', '100', '150', '200','no'})
    xlim([0 270])
    title(stims{stim});
    ylim([400, 1100])
    xlabel('Laser onset (ms)')
    ylabel('Response time (ms)')
end
    

%% do ANOVA

%test for normality
for stim = 1:3
    figure;
    RTtask = meanRTsAll{stim};
    %remove mice without data
    notnans = ~isnan(RTtask);
    countnonans = sum(notnans,2);
    RTtask = RTtask(countnonans>3,:);
    %remove latencies without data
    notnans = ~isnan(RTtask);
    countnonans = sum(notnans,1);
    RTtask = RTtask(:,countnonans > 3);
    
    %for anova
     nmice = size(RTtask,1);
     nlats = size(RTtask,2);
          
 
    %use anovan 1 way rm anova - get stats 
    groupsM = 1:nmice;
    groupsM = repmat(groupsM',1,nlats);
    mousenums = reshape(groupsM,[],1);
    groupsL = 1:nlats;
    groupsL = repmat(groupsL',1,nmice);
    latys = reshape(groupsL',[],1);
    [PS.avnp{stim},PS.avntbl{stim},PS.avnstats{stim}] = anovan(reshape(RTtask,[],1),{latys,mousenums},'random',2);
    PS.multcomp{stim} = multcompare(PS.avnstats{stim});

%     % bonferroni correction
%     PS.OptoP(stim,:) = Pl* (size(RTtask,2)-1);
    
end


    


end
