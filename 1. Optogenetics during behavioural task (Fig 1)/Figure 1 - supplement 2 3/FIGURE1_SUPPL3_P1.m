
% FIGURE1_SUPPL3_P1_November2023 (OptoStim figs)
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
indResponsive = [];
matPsth_OS = [];
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
    indResponsive = [indResponsive sLoad.vecResp(:)'];
    matPsth_OS = [matPsth_OS sLoad.matPsth_OS'];
end
vecTime = sLoad.vecTime; %for PSTHs
vecTimeOS = sLoad.vecTime_OS;

%% GET INDICES
%make sure units come from BL6 mice
indGad2Cre = contains(cellSubjectType,'gad2cre','IgnoreCase',true);
vecBL6RecIdx = find(~indGad2Cre);
indBL6Units = ismember(vecRecIdx,vecBL6RecIdx);
indRem = vecMeanRate_OptoOFF_ALL == 0 & vecMeanRate_OptoON_ALL == 0;
indBL6Units(indRem) = 0;

%%
figure; maxfig;
intRows = 2; intCols = 3;

%% PSTH
subplot(intRows,intCols,2); hold on
vecTimeOS_ms = vecTimeOS*1000; %ms
vecMeanRateShort = mean(matPsth_OS(:,indBL6Units & ~indResponsive)','omitnan'); %#ok<UDIM>
vecSemRateShort = std(matPsth_OS(:,indBL6Units & ~indResponsive)','omitnan')/sqrt(sum(~isnan(matPsth_OS(1,indBL6Units & ~indResponsive)))); %#ok<UDIM>

p1 = shadedErrorBar(vecTimeOS_ms,vecMeanRateShort,vecSemRateShort,'lineProps','k--');
p1.patch.LineStyle = 'none';
p1.edge(1).LineStyle='none';
p1.edge(2).LineStyle=' none';

vecMeanRateShort = mean(matPsth_OS(:,indBL6Units & indResponsive)','omitnan'); %#ok<UDIM> 
vecSemRateShort = std(matPsth_OS(:,indBL6Units & indResponsive)','omitnan')/sqrt(sum(~isnan(matPsth_OS(1,indBL6Units & indResponsive)))); %#ok<UDIM> 
p2 = shadedErrorBar(vecTimeOS_ms,vecMeanRateShort,vecSemRateShort,'lineProps','k');
p2.patch.LineStyle = 'none';
p2.edge(1).LineStyle='none'; p2.edge(2).LineStyle='none';
xline(0,'b-');
xline(50,'b-');
legend([p2.mainLine p1.mainLine],{'Visual','Non-vis.'});
xlim([-50 250])
xticks([0 100 200]);
yticks([0 2 4 6]);
xlabel('Time from laser onset (ms)');
ylabel('Rate (sp/s)');
fixfig;

%% SCATTER PLOTS MEAN RATE DURING LASER VS BASELINE
%option 1: log-log
subplot(intRows,intCols,3); hold on
vecXref = 10.^(-1:.1:3);
vecYref = vecXref;
plot(vecXref,vecYref,'k--');
s1 = scatter((vecMeanRateBase_OS(~indResponsive & indBL6Units)+1),(vecMeanRateStim_OS(~indResponsive & indBL6Units)+1),40,'k','o','jitter','on');
s2 = scatter((vecMeanRateBase_OS(indResponsive & indBL6Units)+1),(vecMeanRateStim_OS(indResponsive & indBL6Units)+1),40,'filled','k','o','MarkerEdgeColor','k','jitter','on');
set(gca,'YScale','log','XScale','log');
xlabel('Rate + 1 laser off -200-0ms (sp/s)');
ylabel('Rate + 1 laser on 0-200ms (sp/s)');
legend([s2 s1],{'Visual','Non-vis.'},'location','northwest')
xlim([10^0 10^2]);
ylim([10^0 10^2]);
axis square;
fixfig;
