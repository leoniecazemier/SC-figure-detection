%% load data and predefine some stuff
load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat')
% load('N:\Ephys\Analysis\behTanks\eyeDat.mat');
load('C:\Users\leoni\Documents\NIN\Data\eyeDat.mat')


%%  INFO/SETTINGS
% trialgroups:
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

px = -0.1995:0.001:1.7495;
pxl = -1.4995:0.001:0.9995;
nu = size(unit,2);
stime = 0.05;
mtime = 0.1;
qtime = 0.2;
et = px>stime & px<= mtime; %early time vis
lt = px>mtime & px<=qtime; % late time vis
pret = pxl>-0.5 & pxl <= 0; %early time lick
post = pxl >0  & pxl <=0.5;  %late time lick
base = px> -0.15 & px<0;
% whichUnits = find([unit.use]);
%ncUse = length(whichUnits);

colors = {[0,0.9,0],[0,0,1],[1,0.6,0],[0.5,0,0.8]};
trtypes = [1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18];  %figori grori figoop groop figgrey grgrey
rfs = {'inside','outside', 'edge', 'unclear'};
fg = {'figure', 'ground'};
allmice = unique({unit.mouse});
miceid = {unit.mouse};
mousecolormap;
smfvis = 80;
dprimewin = find(px>stime & px<=qtime);
mpetime = 0.65;
mpltime = 0.90;
mltime = 1.161;
lplwin = pxl > -0.7 & pxl <= -0.2;
lpvfwin = px>0 & px < mltime;
lpvwin =  px > mpetime & px <= mpltime;
contwin1 = px > 0 & px <= 0.4;
contwin2 = px > 1.1 & px <= 1.5;

allses = unique({unit.ses});
sessions = unique({unit.ses});
sesList = {unit.ses};

nB = 1000;  %number of bootstraps
mpcs = find([unit.MP]);

% Get mean lick time
RTS = {unit.RTs};
RTS = RTS(mpcs);
mRT = nan(1,length(RTS));
for cs = 1:length(RTS)
    mRT(cs) = nanmean(RTS{cs});
end

meanRT = mean(mRT)/1000;

%% pick which analysis type

mysel = 'MP';
myresp = 'nrDataLick'; % or nrDataLick
mysplits = {'aa'}; % 'CI'}; %CI contra/ipsi or FG figure/ground
stims = 1:2; % Ori, OOP, Base, Ori/OOP

stimnames = {'Ori', 'OOP', 'Grey','Ori OOP'};
stimsels = { [unit.useOri],[unit.useOOP],([unit.useOri] | [unit.useOOP]), ([unit.useOri] | [unit.useOOP])};

% colors = {[0.7 0 1],[0.2 0 0.3],[1 0.5 0],[0.5 0.15 0]};
colors = {[235 32 39]./256,[121 10 10]./256,[68 150 200]./256,[31 67 132]./256,[0 210 75]./256,[0 102 51]./256};

trHelp = [0 2 4];

%% plot eye movements
% y position is not recorded well so only use x position

selection = find([unit.use] & ([unit.useOri] | [unit.useOOP]) & [unit.(mysel)]);

if strcmp(myresp, 'nrDataLick')
    mywin = lplwin;
    mypx = pxl;
    atime = mpetime;
    btime = mpltime;
    myxlim = [-1.5 1];
    myylim = [-1.5 2.5];
    myxpos = 'xposL';
    myypos = 'yposL';
elseif strcmp(myresp, 'nrDataVis')
    mywin = lpvfwin;
    mypx = px;
    atime = 0;
    btime = mltime;
    myxlim = [-0.15 1.7];
    myylim = [-1 3];
    myxpos = 'xpos';
    myypos = 'ypos';
end

for uu=1:length(selection)
    uu = selection(uu);
    tg = unit(uu).trialGroupsNE;
    sesno  = find(strcmp(unit(uu).ses, allses));
    figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8]);
    
    %plot things...
    subplot(2,3,[1 2])
    smdataF = squeeze(dataH{4}(1,uu,:));
    smdataG = squeeze(dataH{4}(2,uu,:));
    plot(mypx, smdataF, 'color',colors{1});
    hold on;
    plot(mypx, smdataG,'color', colors{2});
    xlim(myxlim)
    ylim(myylim)
    yLimits = get(gca, 'ylim');
    stimline = line([0 0], [-2 4], 'color','k', 'linewidth',1.6, 'HandleVisibility','off');
    %     pat = patch([mpetime mpltime mpltime mpetime ],[yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.5 0.5 0.5], 'handleVisibility','off');
    %     alpha(pat,0.2);
    %     set(pat,'EdgeColor','none');
    title(['Ori OOP Hit unit ' num2str(uu)]);
    xlabel('Time from lick (s)');
    ylabel('Normalized firing rate');
    legend({'Figure', 'Ground'});
    
    
    subplot(2,3,[4 5])
    smdataF = squeeze(dataE{4}(1,uu,:));
    smdataG = squeeze(dataE{4}(2,uu,:));
    plot(mypx, smdataF, 'color',colors{1});
    hold on;
    plot(mypx, smdataG,'color', colors{2});
    xlim(myxlim)
    ylim(myylim)
    yLimits = get(gca, 'ylim');
    stimline = line([0 0], [-2 4], 'color','k', 'linewidth',1.6, 'HandleVisibility','off');
    %     pat = patch([mpetime mpltime mpltime mpetime ],[yLimits(1) yLimits(1) yLimits(2) yLimits(2)], [0.5 0.5 0.5], 'handleVisibility','off');
    %     alpha(pat,0.2);
    %     set(pat,'EdgeColor','none');
    title(['Ori OOP Error unit ' num2str(uu)]);
    xlabel('Time from lick (s)');
    ylabel('Normalized firing rate');
    legend({'Figure', 'Ground'});
    
    if unit(uu).eyedata
        %x pos hit
        subplot(2,3,3)
        %             contralicks = find(ismember(tg, [1 5]));
        %             ipsilicks = find(ismember(tg, [2 4]));
        Figs = find(ismember(tg, [1 2 3 7 8 9])); %ori and oop together
        grlicksH = find(ismember(tg, [4 5 6 10 11 12]));
        stimline = line([0 0], [-2 4], 'color','k', 'linewidth',1.6, 'HandleVisibility','off');
        hold on;
        xp1 = nanmean(eyeDat(sesno).(myxpos)(Figs,:),1);
        xp2 = nanmean(eyeDat(sesno).(myxpos)(grlicksH,:),1);
        plot(mypx,xp1, 'color', colors{1});
        plot(mypx,xp2 ,'color', colors{2});
        mymin = min([xp1 xp2],[],'all');
        mymax = max([xp1 xp2], [],'all');
        ylim([mymin-0.05 mymax+0.05])
        xlim(myxlim)
        title(['Ori Hit eye x-position ' num2str(uu)]);
        xlabel('Time from lick (s)');
        ylabel('x-position');
        legend({'Figure', 'Ground'});
        
        
        %x pos error
        subplot(2,3,6)
        %             contralicks = find(ismember(tg, [1 5]));
        %             ipsilicks = find(ismember(tg, [2 4]));
        figlicksE = find(ismember(tg, [2]));
        grlicksE = find(ismember(tg, [5]));
        stimline = line([0 0], [-2 4], 'color','k', 'linewidth',1.6, 'HandleVisibility','off');
        hold on;
        xp1 = nanmean(eyeDat(sesno).(myxpos)(figlicksE,:),1);
        xp2 = nanmean(eyeDat(sesno).(myxpos)(grlicksE,:),1);
        plot(mypx,xp1, 'color', colors{1});
        plot(mypx,xp2 ,'color', colors{2});
        mymin = min([xp1 xp2],[],'all');
        mymax = max([xp1 xp2], [],'all');
        ylim([mymin-0.05 mymax+0.05])
        xlim(myxlim)
        title(['Ori Error eye x-position ' num2str(uu)]);
        xlabel('Time from lick (s)');
        ylabel('x-position');
        legend({'Figure', 'Ground'});
        
    end
end



%% compute correlations between eye and neural response

selection = find([unit.use] & ([unit.useOri] | [unit.useOOP]) & [unit.(mysel)]);
corr_ps = nan(length(selection),1); %correlation for each neuron
corr_coefs = nan(length(selection),1); %correlation for each neuron
%dim 2 is type -  visual/lick, dim 3 is movement type: x/y pos

%make a lowpass filter for correcting eye data glitches
sf = 1017.25;
fpass = 10;
% y = lowpass(x,fpass,fs) specifies that x is sampled at a rate of fs hertz. fpass is the passband frequency of the filter in hertz

useMeans = 0;

ucount = 0;
for uu = selection
    ucount= ucount+1;
    if ~unit(uu).eyedata
        continue
    end
    
    sesno  = find(strcmp(unit(uu).ses, allses));
    
    
    %make mean of all responses
    %(not a vector of all data points because these include a lot of zeros
    %the eye data is also a bit flat/bad range so many same values give an
    %artifically high correlation)
    
    %also make vector of all eye x data, visual
    
    eyeXvecPrep = eyeDat(sesno).(myxpos);
    if useMeans
        respVecV = nanmean(unit(uu).nrDataVis,1);
        eyeXvecV = nanmean(eyeXvecPrep,1);
    else
        respVecV = reshape(unit(uu).nrDataVis',1,[]);
        eyeXvecV = reshape(eyeXvecPrep',1,[]);
    end
    
    %lowpass the eye data because of recoring 'blips' when the tracker lost
    %the pupil, usually because whisker moved
    eyeXvecV = lowpass(eyeXvecV(~isnan(eyeXvecV)),fpass,sf);
    respVecV = respVecV(~isnan(eyeXvecV));
    
    %clean it
    myz = zscore(eyeXvecV);
    notValid = abs(myz) > 5;
    eyeXvecV(notValid) = NaN;
    eyeXvecV = fillmissing(eyeXvecV, 'linear');
    
    %smooth it?
    eyeXvecV = smooth(eyeXvecV,100);
    respVecV = smooth(respVecV,100);
    %make speed from position
    eyeXspeed = diff(eyeXvecV);
    pxspeed = px(2:end);
    respVecSpeed = respVecV(2:end);
    
    figure;
    hold on;
    if useMeans
        plot(pxspeed(20:end-20),zscore(respVecSpeed(20:end-20)), 'k')
        plot(pxspeed(20:end-20),zscore(eyeXspeed(20:end-20)), 'r')
        xlabel('Time from stim onset')
    else
        plot(zscore(respVecSpeed), 'k')
        plot(zscore(eyeXspeed), 'r')
        xlabel('Time')
    end
    ylabel('z-score');
    legend({'Normalized firing rate' 'Eye movement speed'})
    title(['Unit ' num2str(uu)]);
    xlim([20 length(eyeXspeed)-20]);
    ylim([-10 10]);
    
    
    
    %correlating response and x pos eye
    [R, P] = corrcoef(respVecSpeed, eyeXspeed);
    corr_coefs(ucount) = R(1,2);
    corr_ps(ucount) = P(1,2);
    
end

%test whether they are different from chance

myCC = corr_coefs(~isnan(corr_coefs));
[hCCl, pCCl] = lillietest(myCC);
if pCCl>0.05
   [hCoef, pCoef,~,CoefStats] =  ttest(myCC);
else
    warning('This non-parametric test was not coded...');
end


%% Plot example eye movement and summed coefficients for paper

%settings, picked neuron 186
uu = 186;
tg = unit(uu).trialGroupsNE;
sesno  = find(strcmp(unit(uu).ses, allses));
Ts{1,1}= find(ismember(tg, [1 7])); %fig Contra
Ts{1,2} = find(ismember(tg, [2 8])); %Fig Ipsi
Ts{2,1}= find(ismember(tg, [5 11])); %Gr contra
Ts{2,2} = find(ismember(tg, [4 10])); %gr ipsi
YY = [-1.5 3.3];
legLabels = {'Neural activity', 'Eye movement speed'};
colorsE = {[68 150 200]./256,[0.2 0.2 0.2]};
titles = {'Fig Contra lick', 'Fig Ipsi lick';'Ground Contra lick','Ground Ipsi lick'};


for FType = 1:2
    afig = figure;
    pn=0;
    for LType = 1:2
        %get info spikes
        pickTrs = Ts{FType,LType};
        SpikesZ = zscore(unit(uu).nrDataVis(pickTrs,:),0,'all');
        figSem{FType,LType} = movmean(sem(SpikesZ,1),smfvis);
        figMean{FType,LType} = movmean(nanmean(SpikesZ,1),smfvis);
        
        %also get eye position
        Xpos = eyeDat(sesno).xpos(pickTrs,:);
        Xspeed = diff(Xpos,[],2);
        XspeedZ= zscore(Xspeed,0,'all');
        XspeedMean{FType,LType} = movmean(nanmean(XspeedZ,1),smfvis);
        XspeedSem{FType,LType} = movmean(sem(XspeedZ,1),smfvis);
        
        %plot it
        %figure trials
        pn = pn+1;
        subplot(1,2,pn)
        hold on;
        [hFfillSpikes(pn), hFlineSpikes(pn)] = errorfill(px(2:end),figMean{FType,LType}(2:end), figSem{FType,LType}(2:end), colorsE{1});
        [hFfillEye(pn), hFlineEye(pn)] = errorfill(px(2:end),XspeedMean{FType,LType}, XspeedSem{FType,LType}, colorsE{2});
        legend([hFlineSpikes(pn) hFlineEye(pn)] ,legLabels)
        xlabel('Time from stim onset (s)');
%         stimline = line([0 0], [-2 4], 'color','k', 'linewidth',1.6,'HandleVisibility','off');
        axis tight
        ylim(YY)
        ylabel('z-score');
        title(titles{FType,pn});
    end
end

%plot coefficients
figure;
hold on;
Coef = corr_coefs(~isnan(corr_coefs));
semCoef = sem(Coef);
bar(1, mean(Coef),'facecolor',colorsE{1})
errorbar(mean(Coef), semCoef, 'k')
ylim([-0.15 0.15]);
ax1 = gca;                  
ax1.XAxis.Visible = 'off'; 
%add individual values
mps = find([unit.MP]);
mpsIncluded = mps(~isnan(corr_coefs));
exIdx = find(mpsIncluded == 186);
Groups = ones(1, length(mpsIncluded));
Groups(exIdx) = 2;
gscatter(ones(1,length(Groups)), Coef, Groups, [[0 0 0];[1 0 0]] ,'.',20, 'off');
ylabel('Correlation coefficient')
