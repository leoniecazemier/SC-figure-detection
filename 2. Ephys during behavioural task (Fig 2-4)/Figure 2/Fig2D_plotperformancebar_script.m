%% 
% cd N:\Ephys\Analysis\LogsEphysBehaviour
cd C:\Users\leoni\Documents\NIN\Data\LogsEphysBehaviour
files = fullfile('*.mat');
logfiles = dir(files);

%load info on good/bad channels per session
% load('N:\Ephys\Analysis\behTanks\recInfo.mat')
load('C:\Users\leoni\Documents\NIN\Data\recInfo.mat')
% load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat')
load('C:\Users\leoni\Documents\NIN\Data\unitdatareviewed.mat')

%order is fighit, figerror, grndhit, grnderror
%prep matrix for 40 depths and 4 trial cases

%%
for T = 1:length(logfiles)
    fullfilename = logfiles(T).name;
    Tank(T).Name = fullfilename(1:end-16);
    Tank(T).Date = fullfilename(end-14:end-7);
    Tank(T).Tankname = [Tank(T).Name '_' Tank(T).Date];
    Tank(T).Blockno = fullfilename(end-4);
    B = Tank(T).Blockno;
    
    filenameload = ['N:\Ephys\Analysis\BehTanks\Tank' Tank(T).Tankname '_B' B 'ChSel.mat'];
    load(filenameload)
   
    figHitOri(T)   = length(myTank.fighitOri);
    figHitOOP(T)   = length(myTank.fighitOOP);
    figErrorOri(T) = length(myTank.figerrorOri);
    figErrorOOP(T) = length(myTank.figerrorOOP);
    grndHitOri(T)  = length(myTank.grndhitOri);
    grndHitOOP(T)  = length(myTank.grndhitOOP);
    grndErrorOri(T) = length(myTank.grnderrorOri);
    grndErrorOOP(T) = length(myTank.grnderrorOOP);
    figHitGrey(T)   = length(myTank.fighitGrey);
    figErrorGrey(T) = length(myTank.figerrorGrey);
    grndHitGrey(T)  = length(myTank.grndhitGrey);
    grndErrorGrey(T) = length(myTank.grnderrorGrey);

    perfOri(T) = (figHitOri(T)+grndHitOri(T))/...
        (figHitOri(T)+ figErrorOri(T)+ grndHitOri(T)+grndErrorOri(T));
    perfOOP(T) = (figHitOOP(T)+grndHitOOP(T))/...
        (figHitOOP(T)+ figErrorOOP(T)+ grndHitOOP(T)+grndErrorOOP(T));   
    perfGrey(T) = (figHitGrey(T)+grndHitGrey(T))/...
        (figHitGrey(T)+ figErrorGrey(T)+ grndHitGrey(T)+grndErrorGrey(T));
    
end

%% fix it a bit
GreySum = figHitGrey + figErrorGrey + grndHitGrey + grndErrorGrey;

noGrey = find(GreySum<=1);
figHitGrey(noGrey) = NaN;
figErrorGrey(noGrey) = NaN;
GrndHitGrey(noGrey) = NaN;
GrndErrorGrey(noGrey) = NaN;

%remove sessions with low performance - not included. 
perfGrey(noGrey) = NaN;

perfGrey = perfGrey.*100;
perfOriDef = perfOri([recInfo.useOri]==1).*100;
perfOOPDef = perfOOP([recInfo.useOOP]==1).*100;


%% Make it per mouse

uniqueMice = unique({Tank.Name});
for m = 1:length(uniqueMice)
    tanksel = find(strcmpi({Tank.Name},uniqueMice(m)));
    %structure is task - stim - resp
    ntrials(m,1,1,1) = sum(figHitGrey(tanksel));
    ntrials(m,1,1,2) = sum(figErrorGrey(tanksel));
    ntrials(m,1,2,1) = sum(grndHitGrey(tanksel));
    ntrials(m,1,2,2) = sum(grndErrorGrey(tanksel));
    ntrials(m,2,1,1) = sum(figHitOri(tanksel));
    ntrials(m,2,1,2) = sum(figErrorOri(tanksel));
    ntrials(m,2,2,1) = sum(grndHitOri(tanksel));
    ntrials(m,2,2,2) = sum(grndErrorOri(tanksel));
    ntrials(m,3,1,1) = sum(figHitOOP(tanksel));
    ntrials(m,3,1,2) = sum(figErrorOOP(tanksel));
    ntrials(m,3,2,1) = sum(grndHitOOP(tanksel));
    ntrials(m,3,2,2) = sum(grndErrorOOP(tanksel));
    
    Performances(m,1) = (ntrials(m,1,1,1)+ntrials(m,1,2,1)) / ...
        (ntrials(m,1,1,1)+ntrials(m,1,2,1)+ ntrials(m,1,1,2)+ntrials(m,1,2,2));
    Performances(m,2) = (ntrials(m,2,1,1)+ntrials(m,2,2,1)) / ...
        (ntrials(m,2,1,1)+ntrials(m,2,2,1)+ ntrials(m,2,1,2)+ntrials(m,2,2,2));
    Performances(m,3) = (ntrials(m,3,1,1)+ntrials(m,3,2,1)) / ...
        (ntrials(m,3,1,1)+ntrials(m,3,2,1)+ ntrials(m,3,1,2)+ntrials(m,3,2,2));
end
    
Performances = Performances*100;


%% plot it

figure; 
hold on;
labels = {'Contrast','Orientation','Phase'};
b = bar(nanmean(Performances, 1));
t = errorbar(1:3,nanmean(Performances,1), nansem_large(Performances,1),'.');
ylim([0 100]);
set(gca, 'xtick', [1:3], 'xticklabel', labels, 'linewidth', 1.2, 'fontsize', 14);
line([0 4],[50 50],'color',[0.5 0.5 0.5],'linestyle', '--','linewidth',1.5);
yticks(0:25:100)

for task=1:3
xvals = rand(5,1)*0.4 + task-0.2;
plot(xvals,Performances(:,task),'o','markersize',10,'color',[0.3 0.3 0.3])
end

b.BarWidth = 0.5; 
b.EdgeColor = 'k';
b.LineWidth = 1.2;
t.LineWidth = 1.8;
b.FaceColor = 'flat';
b.CData(1,:) = [0 0.8 0];  %bacontr is green
b.CData(2,:) = [0.8 0 0];  %ori is red
b.CData(3,:) = [0 0 0.8];  %phase is blue
t.Color = 'k';
set(gcf, 'Position', [200 300 400 600])
ylabel('Accuracy (%)');




