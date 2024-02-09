%% load data and predefine some stuff
load('C:\Users\leoni\Documents\NIN\Data\unitdatareviewed.mat')
% load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat')
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
mpcs = find([unit.MP]);

% Get mean lick time
RTS = {unit.RTs};
RTS = RTS(mpcs);
mRT = nan(1,length(RTS));
for cs = 1:length(RTS)
    mRT(cs) = nanmean(RTS{cs});
end

meanRT = mean(mRT)/1000;

stimNames = {'Orientation','Phase'};
responseNames = {'Hit','Error'};
    


%% pick which analysis type

mysel = 'MP';
myresp = 'nrDataLick'; % or nrDataLick
mysplits = {'aa'}; % 'CI'}; %CI contra/ipsi or FG figure/ground
stims = 1:2; % Ori, OOP, Base, Ori/OOP

stimsels = { [unit.useOri],[unit.useOOP],([unit.useOri] | [unit.useOOP]), ([unit.useOri] | [unit.useOOP])};

% colors = {[0.7 0 1],[0.2 0 0.3],[1 0.5 0],[0.5 0.15 0]};
colors = {[235 32 39]./256,[121 10 10]./256;[68 150 200]./256,[31 67 132]./256};

trHelp = [0 2 4];

%% plot what cells with motor peak are doing

%also plot lick centered firing rates of motor peak neurons


for spl = 1:length(mysplits)
    mysplit = mysplits{spl};
    if strcmpi(mysplit, 'FG')
        trH = {[1 2 3; 4 5 6],[7 8 9;10 11 12],[13 14 15; 16 17 18], [1 2 3 7 8 9;4 5 6 10 11 12]};
        figns = [40 41];
        ttypes = 1:2;
        leg = {'Figure','Ground'};
    elseif strcmpi(mysplit,'CI')
        trH = {[1 5; 2 4],[7 11;8 10],[13 17; 14 16], [1 5 7 11;2 4 8 10]};
        figns = [50 51];
        ttypes = 1:2;
        leg = {'Contra lick','Ipsi lick'};
    else  % all 4 options separate FH FE GH GE
        trH = {[1,4],[2,5];[7,10],[8,11]};
        fgs = 1:2; 
        reactions = 1:2;
        figns = [60 61];
        leg = {'Figure Contra lick (Hit)', 'Figure Ipsi lick (Error)', ...
            'Ground Ipsi lick (Hit)', 'Ground Contra lick (Error)'};
    end
    
    figure(figns(2));
    clf
    plotno=0;
    clear data
    for stim = stims % stim3 is base task but not enough data


        clear smdatab
        clear earlydat
        clear latedat
        clear mymice
        
        hold on
        selection = find([unit.use] & stimsels{stim} & [unit.(mysel)]);
        n = length(selection);
        for reaction = reactions
            plotno=plotno+1;
            subplot(2,2,plotno)
            for fg = fgs
                for c=1:length(selection)
                    uu = selection(c);  % uu is the unit of choice
                    tSelect = find(ismember(unit(uu).trialGroups, trH{stim,reaction}(fg)));
                    if isempty(tSelect)
                        data{stim,reaction,fg}(c,:) = nan(1,2500);
                        nn(ttype) = nn(ttype)-1;
                        fprintf('Skipped cell %i in ori lick plot. \n', uu)
                    else
                        data{stim,reaction,fg}(c,:) = nanmean(unit(uu).nrDataLick(tSelect,:),1);
                    end
                    mymice(c) = find(strcmp(allmice,unit(uu).mouse));
                end
                smdatab(:,ttype) = movmean(nanmean(data{stim,reaction,fg},1),smfvis);
                col = colors{stim,fg};
                semb = movmean(nansem_large(data{stim,reaction,fg},1)',smfvis);
                [handleFill,handleLine] = errorfill(pxl,smdatab(:,ttype),semb,col);
            end
%             stimline = line([0 0],[-2 2] , 'color','k', 'linewidth',1.6, 'HandleVisibility','off');
            xlim([-1 0.5])
            ylim([-1 3]);
            yLimits = get(gca, 'ylim');
            if strcmpi(mysplit,'FG')
                title(['FGM lick centered ' stimnames{stim} ', n=' num2str(min(nn))]);
            elseif strcmpi(mysplit,'CI')
                title(['Contra-ipsi lick - LC ' stimnames{stim} ', n=' num2str(min(nn))]);
            else
                title(['Lick-centered ' stimNames{stim} ', ' responseNames{reaction}])
            end
            
        end
        xlabel('Time(s)');
        ylabel('Normalized firing rate');
        set(gca, 'tickdir', 'out');
        legend({'Figure','Ground'});
        
        
    end
end
