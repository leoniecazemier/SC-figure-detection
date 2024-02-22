% FIGURE1_SUPPL3_P2_November2023 (OptoGratings figs)
%   Supplementary analysis for 'Involvement of superior colliculus in complex figure detection of mice', Cazemier et al.

clear; close all;
runHeaderGAD2;

%% LOAD DATA
cellSubjectType = [];
vecMeanRate_OptoON_ALL = [];
vecMeanRate_OptoOFF_ALL = [];
vecRateSpontaneous_ALL = [];
matPsthON_ALL = [];
matPsthOFF_ALL = [];
vecDepth = [];
vecRecIdx = [];
vecMouseNum = [];
indResponsive = [];
vecZetaP_OS = [];
vecCluIdx = [];
vecMeanRateBase_OS = [];
vecMeanRateStim_OS = [];
if ~exist(strSavePath,'dir'), error('Results path does not exist!'); end

sFiles = dir(fullpath(strSavePath,'OptoGratings*.mat'));
for intFile=1:numel(sFiles)
    fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
    sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
    cellSubjectType{intFile} = sLoad.sJson.subjecttype; %#ok<SAGROW>
    vecMeanRate_OptoON_ALL = [vecMeanRate_OptoON_ALL sLoad.vecMeanRate_OptoON(:)'];
    vecMeanRate_OptoOFF_ALL = [vecMeanRate_OptoOFF_ALL sLoad.vecMeanRate_OptoOFF(:)'];
    vecRateSpontaneous_ALL = [vecRateSpontaneous_ALL sLoad.vecRateSpontaneous(:)'];
    matPsthON_ALL = [matPsthON_ALL sLoad.matPsthON'];
    matPsthOFF_ALL = [matPsthOFF_ALL sLoad.matPsthOFF'];
    vecDepth = [vecDepth sLoad.vecDepth(:)'];
    vecRecIdx = [vecRecIdx repmat(intFile,size(sLoad.vecMeanRate_OptoON(:)'))];
    vecMouseNum = [vecMouseNum repmat(str2double(sLoad.sJson.subject),size(sLoad.vecMeanRate_OptoON(:)'))];
    indResponsive = [indResponsive sLoad.vecResp(:)'];
    vecMeanRateBase_OS = [vecMeanRateBase_OS sLoad.vecMeanRateBase_OS(:)'];
    vecMeanRateStim_OS = [vecMeanRateStim_OS sLoad.vecMeanRateStim_OS(:)'];
    vecZetaP_OS = [vecZetaP_OS sLoad.vecZetaP_OS(:)'];
    vecCluIdx = [vecCluIdx sLoad.vecClusterIdx(:)'];
end
vecTime = sLoad.vecTime; %for PSTHs, same for all recs
indResponsive = logical(indResponsive);

%% GET INDICES
%make sure units come from BL6 mice
indGad2Cre = contains(cellSubjectType,'gad2cre','IgnoreCase',true);
vecBL6RecIdx = find(~indGad2Cre);
indBL6Units = ismember(vecRecIdx,vecBL6RecIdx);
indRem = vecMeanRate_OptoOFF_ALL == 0 & vecMeanRate_OptoON_ALL == 0;
indBL6Units(indRem) = 0;

%% BL6 - POPULATION RESPONSE GRATINGS OPTO ON VS OFF
figure; maxfig;
intRows = 2; intCols = 4;

%PSTH
%include responses from ALL units
vecMeanPSTH_OptoOFF = mean(matPsthOFF_ALL(:,indBL6Units),2);
vecSemPSTH_OptoOFF = std(matPsthOFF_ALL(:,indBL6Units ),[],2)/...
    sqrt(size(matPsthOFF_ALL(:,indBL6Units),2));
vecMeanPSTH_OptoON = mean(matPsthON_ALL(:,indBL6Units),2);
vecSemPSTH_OptoON = std(matPsthON_ALL(:,indBL6Units),[],2)/...
    sqrt(size(matPsthON_ALL(:,indBL6Units),2));

subplot(intRows,intCols,[1 2]); hold on
vecTimeMs = vecTime * 1000; %ms
p1 = shadedErrorBar(vecTimeMs,vecMeanPSTH_OptoOFF,vecSemPSTH_OptoOFF,'lineProps','k');
p1.patch.LineStyle = 'none';
p1.edge(1).LineStyle='none'; p1.edge(2).LineStyle='none';
p2 = shadedErrorBar(vecTimeMs,vecMeanPSTH_OptoON,vecSemPSTH_OptoON,'lineProps','b');
p2.mainLine.Color = lines(1);
p2.patch.FaceColor = lines(1);
p2.patch.LineStyle = 'none';
p2.edge(1).LineStyle='none'; p2.edge(2).LineStyle='none';
xlim([0 0.250]*1000);
legend('Laser off', 'Laser on','Box', 'off')
xlabel('Time from stimulus onset (ms)');
ylabel('Rate (sp/s)');
fixfig;

%scatter
subplot(intRows,intCols,3); hold on
vecXref = 10.^(-1:.1:3);
vecYref = vecXref;
plot(vecXref,vecYref,'k--');
s1 = scatter((vecMeanRate_OptoOFF_ALL(~indResponsive & indBL6Units) +1),(vecMeanRate_OptoON_ALL(~indResponsive & indBL6Units)+1),40,'k','o','jitter','on');
s2 = scatter((vecMeanRate_OptoOFF_ALL(indResponsive & indBL6Units)+1),(vecMeanRate_OptoON_ALL(indResponsive & indBL6Units)+1),40,'filled','k','o','MarkerEdgeColor','k','jitter','on');
intNumPoints = numel(s1.XData)+numel(s2.XData);
set(gca,'YScale','log','XScale','log');
xlabel('Rate + 1 laser off (sp/s)');
ylabel('Rate + 1 laser on (sp/s)');
legend([s2 s1],{'Visual','Non-visual'},'location','northwest')
xlim([10^0 10^2]);
ylim([10^0 10^2]);
axis square
fixfig

%stats, LME
subplot(intRows,intCols,4); hold on
vecMeanRate = [vecMeanRate_OptoOFF_ALL(indBL6Units) vecMeanRate_OptoON_ALL(indBL6Units)];
vecOptoOn = [zeros(size(vecMeanRate_OptoOFF_ALL(indBL6Units))) ones(size(vecMeanRate_OptoON_ALL(indBL6Units)))];
vecSession = [vecRecIdx(indBL6Units) vecRecIdx(indBL6Units)];
vecUnit = [1:numel(vecMeanRate_OptoOFF_ALL(indBL6Units)) 1:numel(vecMeanRate_OptoON_ALL(indBL6Units))];
tblBL6 = table(vecMeanRate(:),vecOptoOn(:),vecSession(:),vecUnit(:),...
    'VariableNames',{'rate','opto','session','unit'});
lme_bl6 = fitlme(tblBL6,'rate~opto+(1|unit)+(opto|session)',...
    'FitMethod','REML','CheckHessian',1,'StartMethod','random');
dblPval_Bl6 = lme_bl6.Coefficients.pValue(2);

%plot rates on vs off
c = lines(1);
scatter(ones(size(vecMeanRate_OptoOFF_ALL(indBL6Units))),vecMeanRate_OptoOFF_ALL(indBL6Units),'jitter','on','MarkerFaceColor','k','MarkerEdgeColor','none');
scatter(2*ones(size(vecMeanRate_OptoON_ALL(indBL6Units))),vecMeanRate_OptoON_ALL(indBL6Units),'jitter','on','MarkerFaceColor',c,'MarkerEdgeColor','none')
scatter(1,mean(vecMeanRate_OptoOFF_ALL(indBL6Units)),3000,'r','_')
scatter(2,mean(vecMeanRate_OptoON_ALL(indBL6Units)),3000,'r','_')
sigstar([1 2],dblPval_Bl6);
xlim([0.5 2.5])
xticks([1 2]);
yticks([0 20 40 60 80]);
xticklabels({'Laser off','Laser on'});
ylabel('Rate (sp/s)');
fixfig;
