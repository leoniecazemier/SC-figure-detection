% FIGURE1_SUPPL2_P1_November2023 (OptoStim figs)
%   Supplementary analysis for 'Involvement of superior colliculus in complex figure detection of mice', Cazemier et al.

clear; close all;
runHeaderGAD2;

%% LOAD DATA
cellSubjectType = [];
vecRecIdx = [];
vecMouseNum = [];
vecCluIdx = [];
vecDepth = [];
vecMeanRateBase_OS = [];
vecMeanRateStim_OS = [];
vecZetaP_OS = [];
vecTtestP_OS = [];
vecSpPerPulse_OS = [];
vecFirstSpLat_OS = [];
vecFirstSpJit_OS = [];
vecMeanRate_OptoON_ALL = [];
vecMeanRate_OptoOFF_ALL = [];

% indResponsive = [];
if ~exist(strSavePath,'dir'), error('Results path does not exist!'); end
sFiles = dir(fullpath(strSavePath,'OptoGratings*.mat'));
for intFile=1:numel(sFiles)
    fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
    sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
    cellSubjectType{intFile} = sLoad.sJson.subjecttype; %#ok<SAGROW>
    vecRecIdx = [vecRecIdx repmat(intFile,size(sLoad.vecMeanRate_OptoON(:)'))];
    vecMouseNum = [vecMouseNum repmat(str2double(sLoad.sJson.subject),size(sLoad.vecMeanRate_OptoON(:)'))];
    vecCluIdx = [vecCluIdx sLoad.vecClusterIdx(:)'];
    vecDepth = [vecDepth sLoad.vecDepth(:)'];
    vecMeanRateBase_OS = [vecMeanRateBase_OS sLoad.vecMeanRateBase_OS(:)'];
    vecMeanRateStim_OS = [vecMeanRateStim_OS sLoad.vecMeanRateStim_OS(:)'];
    vecMeanRate_OptoON_ALL = [vecMeanRate_OptoON_ALL sLoad.vecMeanRate_OptoON(:)'];
    vecMeanRate_OptoOFF_ALL = [vecMeanRate_OptoOFF_ALL sLoad.vecMeanRate_OptoOFF(:)'];
    vecZetaP_OS = [vecZetaP_OS sLoad.vecZetaP_OS(:)'];
    vecTtestP_OS = [vecTtestP_OS sLoad.vecttestP_OS(:)'];
    vecSpPerPulse_OS = [vecSpPerPulse_OS sLoad.vecSpPerPulseALL(:)'];
    vecFirstSpLat_OS = [vecFirstSpLat_OS sLoad.vecFirstSpLatALL(:)'];
    vecFirstSpJit_OS = [vecFirstSpJit_OS sLoad.vecFirstSpJitALL(:)'];
%     indResponsive = [indResponsive sLoad.vecResp(:)'];
end

%% GET PUTATIVE(!) GAD2+ UNITS
%make sure units come from GAD2-Cre mice (extra check) 
indGad2Cre = contains(cellSubjectType,'gad2cre','IgnoreCase',true);
vecGad2CreRecIdx = find(indGad2Cre);
indGad2CreUnits = ismember(vecRecIdx,vecGad2CreRecIdx);

indRem = vecMeanRate_OptoOFF_ALL == 0 & vecMeanRate_OptoON_ALL == 0;
indGad2CreUnits(indRem) = 0;
%criteria for 'tagging', based on ±10ms around laser onset:
indPutGad2 = vecMeanRateStim_OS > vecMeanRateBase_OS & vecZetaP_OS < 0.05 & indGad2CreUnits; 
% indPutGad2T = vecMeanRateStim_OS > vecMeanRateBase_OS &vecTtestP_OS < 0.01;

%also get units resp. to laser but not tagged (i.e., downreg.)
indResp = vecZetaP_OS < 0.05 & ~indPutGad2 & indGad2CreUnits;
indUnResp = vecZetaP_OS >= 0.05 & indGad2CreUnits;

%%
figure; maxfig;
intRows = 2; intCols = 3;

%% SCATTER PLOTS MEAN RATE DURING LASER VS BASELINE
c = lines(2);

%option 1: log-log
subplot(intRows,intCols,3); hold on
vecXref = 10.^(-1:.1:3);
vecYref = vecXref;
plot(vecXref,vecYref,'k--');
indExamp =  (vecRecIdx == 8 & vecCluIdx == 435) | (vecRecIdx == 3 & vecCluIdx == 562);
s1 = scatter((vecMeanRateBase_OS(indUnResp)+1),(vecMeanRateStim_OS(indUnResp)+1),40,'k','o','jitter','on','jitteramount',0.1);
s2 = scatter((vecMeanRateBase_OS(indResp)+1),(vecMeanRateStim_OS(indResp)+1),40,'filled','k','o','MarkerEdgeColor','k','jitter','on','jitteramount',0.1);
s3 = scatter((vecMeanRateBase_OS(indPutGad2 & ~indExamp)+1),(vecMeanRateStim_OS(indPutGad2 & ~indExamp)+1),40,'filled','MarkerEdgeColor',c(2,:),'MarkerFaceColor',c(2,:),'jitter','on','jitteramount',0.1);
%example units
indEx1 = vecRecIdx == 8 & vecCluIdx == 435;
indEx2 = vecRecIdx == 3 & vecCluIdx == 562;
scatter((vecMeanRateBase_OS(indEx1)+1),(vecMeanRateStim_OS(indEx1)+1),40,'filled','MarkerFaceColor','c','MarkerEdgeColor','c');
scatter((vecMeanRateBase_OS(indEx2)+1),(vecMeanRateStim_OS(indEx2)+1),40,'filled','MarkerFaceColor','m','MarkerEdgeColor','m');
set(gca,'YScale','log','XScale','log');
xlabel('Rate + 1 laser off -10-0ms (sp/s)');
ylabel('Rate + 1 laser on 0-10ms (sp/s)');
legend([s1 s2 s3],{'p < 0.05','p > 0.05','put. GAD2+'},'location','best')
xlim([10^0 inf]);
ylim([10^0 inf])
axis square
fixfig

% %option 2: only y-log
% subplot(intRows,intCols,3); hold on
% vecXref = 10.^(-1:.1:3);
% vecYref = vecXref;
% plot(vecXref,vecYref,'k--');
% s1 = scatter(vecMeanRateBase_OS(indUnResp),(vecMeanRateStim_OS(indUnResp)+1),40,'k','o','jitter','on');
% s2 = scatter(vecMeanRateBase_OS(indResp),(vecMeanRateStim_OS(indResp)+1),40,'filled','k','o','MarkerEdgeColor','k','jitter','on');
% s3 = scatter(vecMeanRateBase_OS(indPutGad2),(vecMeanRateStim_OS(indPutGad2)+1),40,'filled','MarkerEdgeColor',c,'MarkerFaceColor',c,'jitter','on');
% %example units
% indEx1 = vecRecIdx == 8 & vecCluIdx == 435;
% indEx2 = vecRecIdx == 3 & vecCluIdx == 562;
% scatter((vecMeanRateBase_OS(indEx1)+1),(vecMeanRateStim_OS(indEx1)+1),40,'filled','MarkerFaceColor','c','MarkerEdgeColor','c');
% scatter((vecMeanRateBase_OS(indEx2)+1),(vecMeanRateStim_OS(indEx2)+1),40,'filled','MarkerFaceColor','m','MarkerEdgeColor','m');
% set(gca,'YScale','log');
% xlabel('Rate laser off (sp/s)');
% ylabel('Rate + 1 laser on (sp/s)');
% legend([s1 s2 s3],{'p < 0.05','p > 0.05','put. GAD2+'},'location','best')
% xlim([0 15]);
% ylim([10^0 inf])
% axis square
% fixfig

%% PLOT SOME METRICS OF PUTATIVE GAD2+ neurons
subplot(intRows,intCols,4); hold on %spikes/pulse
vecBinEdges = 0:0.5:6;
histogram(vecSpPerPulse_OS(indPutGad2),vecBinEdges);
% histogram(vecSpPerPulse(indPutGad2T),vecBinEdges);
dblThisMean = mean(vecSpPerPulse_OS(indPutGad2),'omitnan');
xline(dblThisMean,'r');
text(dblThisMean+(max(vecBinEdges)/50),35,num2str(round(dblThisMean,2)),...
    'FontSize',14,'Color','r')
ylabel('Number of neurons');
xlabel('Spikes evoked per pulse');
ylim([0 40]);
yticks([0 10 20 30 40]);
% legend({'zeta','ttest'})
axis square;
fixfig;

subplot(intRows,intCols,5); hold on %first spike latency
vecBinEdges = 0:1:10;
histogram(vecFirstSpLat_OS(indPutGad2)*1000,vecBinEdges);
% histogram(vecFirstSpLat(indPutGad2T)*1000,vecBinEdges);
dblThisMean = mean(vecFirstSpLat_OS(indPutGad2)*1000,'omitnan');
xline(dblThisMean,'r');
text(dblThisMean+(max(vecBinEdges)/50),35,num2str(round(dblThisMean,2)),...
    'FontSize',14,'Color','r')
ylabel('Number of neurons');
xlabel('First spike latency (ms)');
ylim([0 40]);
yticks([0 10 20 30 40]);
axis square;
fixfig;

subplot(intRows,intCols,6); hold on %first spike jitter
vecBinEdges = 0:0.5:5;
histogram(vecFirstSpJit_OS(indPutGad2)*1000,vecBinEdges);
% histogram(vecFirstSpJit(indPutGad2T)*1000,vecBinEdges);
dblThisMean = mean(vecFirstSpJit_OS(indPutGad2)*1000,'omitnan');
xline(dblThisMean,'r');
text(dblThisMean+(max(vecBinEdges)/50),35,num2str(round(dblThisMean,2)),...
    'FontSize',14,'Color','r')
ylabel('Number of neurons');
xlabel('First spike jitter (ms)');
ylim([0 40]);
yticks([0 10 20 30 40]);
axis square;
fixfig;
