% FIGURE1_SUPPL2_P2_November2023 (OptoGratings figs)
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

c = lines(1);
%% GET INDICES
%make sure units come from GAD2-Cre mice (extra check)
indGad2Cre = contains(cellSubjectType,'gad2cre','IgnoreCase',true);
vecGad2CreRecIdx = find(indGad2Cre); indGad2CreUnits = ismember(vecRecIdx,vecGad2CreRecIdx);

%criteria for 'tagging', based on ±10ms around laser onset:
indPutGad2 = vecMeanRateStim_OS > vecMeanRateBase_OS & vecZetaP_OS < 0.05 & indGad2CreUnits;
% indPutGad2T = vecMeanRateStim_OS > vecMeanRateBase_OS &vecTtestP_OS < 0.01;

indRem = vecMeanRate_OptoOFF_ALL == 0 & vecMeanRate_OptoON_ALL == 0;

indPutGad2(indRem) = 0;
indGad2CreUnits(indRem) = 0;
indUnResp = ~indResponsive;

indResponsive(indRem) = 0;
indUnResp(indRem) = 0;
%% GAD2Cre - POPULATION RESPONSE GRATINGS OPTO ON VS OFF
figure;
maxfig;

fprintf('%d units from Gad2Cre animals!\n',sum(indGad2CreUnits));

intRows = 2; intCols = 4;

%psth
%only responses from visually responsive units!
% vecMeanPSTH_OptoOFF = mean(matPsthOFF_ALL(:,indGad2CreUnits & indResponsive),2);
% vecSemPSTH_OptoOFF = std(matPsthOFF_ALL(:,indGad2CreUnits & indResponsive),[],2)/...
%     sqrt(size(matPsthOFF_ALL(:,indGad2CreUnits & indResponsive),2));
% vecMeanPSTH_OptoON = mean(matPsthON_ALL(:,indGad2CreUnits & indResponsive),2);
% vecSemPSTH_OptoON = std(matPsthON_ALL(:,indGad2CreUnits & indResponsive),[],2)/...
%     sqrt(size(matPsthON_ALL(:,indGad2CreUnits & indResponsive),2));

%include responses from ALL units
vecMeanPSTH_OptoOFF = mean(matPsthOFF_ALL(:,indGad2CreUnits),2);
vecSemPSTH_OptoOFF = std(matPsthOFF_ALL(:,indGad2CreUnits ),[],2)/...
    sqrt(size(matPsthOFF_ALL(:,indGad2CreUnits),2));
vecMeanPSTH_OptoON = mean(matPsthON_ALL(:,indGad2CreUnits),2);
vecSemPSTH_OptoON = std(matPsthON_ALL(:,indGad2CreUnits),[],2)/...
    sqrt(size(matPsthON_ALL(:,indGad2CreUnits),2));

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

vecMeanEvokedRate = [vecMeanRate_OptoOFF_ALL(indGad2CreUnits) vecMeanRate_OptoON_ALL(indGad2CreUnits)];
vecOptoOn = [zeros(size(vecMeanRate_OptoOFF_ALL(indGad2CreUnits))) ones(size(vecMeanRate_OptoON_ALL(indGad2CreUnits)))];
vecSession = [vecRecIdx(indGad2CreUnits) vecRecIdx(indGad2CreUnits)];
vecUnit = [1:numel(vecMeanRate_OptoOFF_ALL(indGad2CreUnits)) 1:numel(vecMeanRate_OptoON_ALL(indGad2CreUnits))];
tblGad2 = table(vecMeanEvokedRate(:),vecOptoOn(:),vecSession(:),vecUnit(:),...
    'VariableNames',{'rate','opto','session','unit'});
lme_gad2 = fitlme(tblGad2,'rate~opto+(1|unit)+(opto|session)',...
    'FitMethod','REML','CheckHessian',1,'StartMethod','random');
dblPval_Gad2 = lme_gad2.Coefficients.pValue(2);

%add inset
axes('Position',[.35 .65 .05 .2]); hold on
box on;
c = lines(1);
scatter(ones(size(vecMeanRate_OptoOFF_ALL(indGad2CreUnits))),vecMeanRate_OptoOFF_ALL(indGad2CreUnits),'jitter','on','MarkerFaceColor','k','MarkerEdgeColor','none');
scatter(2*ones(size(vecMeanRate_OptoON_ALL(indGad2CreUnits))),vecMeanRate_OptoON_ALL(indGad2CreUnits),'jitter','on','MarkerFaceColor',c,'MarkerEdgeColor','none')
scatter(1,mean(vecMeanRate_OptoOFF_ALL(indGad2CreUnits)),800,'r','_')
scatter(2,mean(vecMeanRate_OptoON_ALL(indGad2CreUnits)),800,'r','_')
sigstar([1 2],dblPval_Gad2);
xlim([0.5 2.5])
xticks([]);
ylabel('Rate (sp/s)')

%scatter
subplot(intRows,intCols,3); hold on
vecXref = 10.^(-1:.1:3);
vecYref = vecXref;
plot(vecXref,vecYref,'k--');
c = lines(2);
s1 = scatter((vecMeanRate_OptoOFF_ALL(indUnResp & indGad2CreUnits & ~indPutGad2)+1),(vecMeanRate_OptoON_ALL(indUnResp & indGad2CreUnits & ~indPutGad2)+1),40,'k','o','jitter','on');
s2 = scatter((vecMeanRate_OptoOFF_ALL(indResponsive & indGad2CreUnits & ~indPutGad2)+1),(vecMeanRate_OptoON_ALL(indResponsive & indGad2CreUnits & ~indPutGad2)+1),40,'filled','k','o','MarkerEdgeColor','k','jitter','on');
s3 = scatter((vecMeanRate_OptoOFF_ALL(indUnResp & indPutGad2)+1),(vecMeanRate_OptoON_ALL(indUnResp & indPutGad2)+1),40,'MarkerEdgeColor',c(2,:),'jitter','on');
s4 = scatter((vecMeanRate_OptoOFF_ALL(indResponsive & indPutGad2)+1),(vecMeanRate_OptoON_ALL(indResponsive & indPutGad2)+1),40,'filled','MarkerEdgeColor',c(2,:),'MarkerFaceColor',c(2,:),'jitter','on');
intNumPoints = numel(s1.XData)+numel(s2.XData)+numel(s3.XData)+numel(s4.XData);

%example units
indEx1 = vecRecIdx == 8 & vecCluIdx == 435;
indEx2 = vecRecIdx == 3 & vecCluIdx == 562;
scatter((vecMeanRate_OptoOFF_ALL(indEx1)+1),(vecMeanRate_OptoON_ALL(indEx1)+1),40,'filled','MarkerFaceColor','c','MarkerEdgeColor','c');
scatter((vecMeanRate_OptoOFF_ALL(indEx2)+1),(vecMeanRate_OptoON_ALL(indEx2)+1),40,'filled','MarkerFaceColor','m','MarkerEdgeColor','m');
set(gca,'YScale','log','XScale','log');
xlabel('Rate + 1 laser off (sp/s)');
ylabel('Rate + 1 laser on (sp/s)');
legend([s2 s1 s4 s3],{'Visual.','Non-visual','Put. Gad2+ & vis.','Put. Gad2+ & non-vis.'},'location','best')
xlim([10^0 10^2]);
ylim([10^0 10^2]);
axis square
fixfig
%%
%pct change in rate
vecPctChangeRate = -((vecMeanRate_OptoOFF_ALL-vecMeanRate_OptoON_ALL)./vecMeanRate_OptoOFF_ALL)*100;
vecBinEdges = -100:10:100;
vecBinCents = vecBinEdges(1:end-1)+mean(diff(vecBinEdges))/2;
matCounts = nan(4,numel(vecBinEdges)-1);
matCounts(1,:) = histcounts(vecPctChangeRate(indUnResp & indGad2CreUnits & ~indPutGad2),vecBinEdges);
matCounts(2,:) = histcounts(vecPctChangeRate(indResponsive & indGad2CreUnits & ~indPutGad2),vecBinEdges);
matCounts(3,:) = histcounts(vecPctChangeRate(indUnResp & indGad2CreUnits & indPutGad2),vecBinEdges);
matCounts(4,:) = histcounts(vecPctChangeRate(indResponsive & indGad2CreUnits & indPutGad2),vecBinEdges);

%add outlier
vecBinCents = [vecBinCents 105];
vecExtrCount = [1 0 0 0];
matCounts = [matCounts vecExtrCount(:)];

subplot(intRows,intCols,4); hold on
yline(0,'k--','Alpha',1)
b = barh(vecBinCents,matCounts','stacked','BarWidth',0.75);
xline(0)
c = lines(2);
b(1).FaceColor = 'w'; b(1).EdgeColor = 'k';
b(2).FaceColor = 'k'; b(2).EdgeColor = 'k';
b(3).FaceColor = 'w'; b(3).EdgeColor = c(2,:);
b(4).FaceColor = c(2,:); b(4).EdgeColor = c(2,:);
yticks([-100 -50 0 50 100 105])
yticklabels({'-100','-50','0','50','100','> 100'})
xlabel('Number of neurons');
ylabel('Δ Rate (%)');
ylim([-100 110]);
fixfig;

fprintf('Of the units from Gad2Cre animals, %d were upregulated and %d were downregulated, %d did not change\n',...
    sum(vecPctChangeRate(indGad2CreUnits)>0),sum(vecPctChangeRate(indGad2CreUnits)<0),sum(vecPctChangeRate(indGad2CreUnits)==0))

%% Plot delta rate vs sessions / depth
vecThisDepths = vecDepth(indGad2CreUnits);
vecThisRecs = vecRecIdx(indGad2CreUnits);
vecThisUnRecs = unique(vecThisRecs);
vecThisSubs = vecMouseNum(indGad2CreUnits);
vecThisPctChange = vecPctChangeRate(indGad2CreUnits);

%normalize depth to that of the most superficial SC neuron in the rec
vecNormDepths = [];
for intRec = 1:numel(vecThisUnRecs)
    intThisRec = vecThisUnRecs(intRec);
    indThis = vecThisRecs == intThisRec;
    dblThisDepth = (vecThisDepths(indThis) - min(vecThisDepths(indThis)))*1000; %positive
    vecNormDepths = [vecNormDepths dblThisDepth];
end

%plot dRate for 79154 & 79155 w/multiple recs
vecSubsID = [79154 79155];
% c = lines(2);

subplot(intRows,intCols,[5 6]); hold on

for intSub = 1:2
    vecThisThisRecs = ...
        vecThisRecs(vecThisSubs == vecSubsID(intSub));
    vecThisThisChange = ...
        vecThisPctChange(vecThisSubs == vecSubsID(intSub));
    vecThisThisChange(vecThisThisChange>100) = 110;
    vecThisThisUnRec = unique(vecThisThisRecs);
    vecMeanChange = [];
    vecMedianChange = [];
    for intRec = 1:numel(vecThisThisUnRec)
        vecMeanChange(intRec) = mean(vecThisThisChange(vecThisThisRecs == vecThisThisUnRec(intRec))); %#ok<SAGROW>
        vecMedianChange(intRec) = median(vecThisThisChange(vecThisThisRecs == vecThisThisUnRec(intRec))); %#ok<SAGROW> 
    end
    if intSub == 1
        s1 = scatter(vecThisThisRecs,vecThisThisChange,40,...
            'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
            'Jitter','on','JitterAmount',0.1);
        scatter(vecThisThisUnRec,vecMeanChange,500,'r','_')
        scatter(vecThisThisUnRec,vecMedianChange,500,'k','_')
    else
        s2 = scatter(vecThisThisRecs,vecThisThisChange,40,...
            'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
            'Jitter','on','JitterAmount',0.1);
        scatter(vecThisThisUnRec,vecMeanChange,500,'r','_')
        scatter(vecThisThisUnRec,vecMedianChange,500,'k','_')

    end
end

% legend([s1 s2],{'Mouse 79154','Mouse 79155'},'Location','best')
xticks(1:1:8);
xlim([0.5 8.5])
xticklabels({'1','2','3','4','1','2','3','4'});
yticks([-100 -50 0 50 100 110]);
yticklabels({'-100','-50','0','50','100','> 100'})
ylim([-100 110])
xlabel('Session');
ylabel('Δ Rate (%)');
fixfig;

%

vecMeanChange2 = [];
vecMedianChange2 = [];
vecEdges = [0 100 200 300 400 500];
cellRateChange = {};
for intBin = 1:numel(vecEdges)-1
    indThis = vecNormDepths >= vecEdges(intBin) & vecNormDepths < vecEdges(intBin + 1);
    dblMedianChange = median(vecThisPctChange(indThis));
    dblMeanChange = mean(vecThisPctChange(indThis));
    cellRateChange{intBin} = vecThisPctChange(indThis); %#ok<SAGROW>
    vecMeanChange2 = [vecMeanChange2 dblMeanChange];
    vecMedianChange2 = [vecMedianChange2 dblMedianChange];
end

subplot(intRows,intCols,[7.5 8]); hold on
yline(-100,'k--','Alpha',.5)
vecM = vecEdges(1:end-1)+50; % for plotting
for intBin = 1:numel(vecM)
    r = cellRateChange{intBin};
    r(r>100) = 110;
    scatter(vecM(intBin)*ones(size(r)),r,40,...
        'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none',...
        'Jitter','on','JitterAmount',10);
end
scatter(vecM,vecMeanChange2,500,'r','_')
scatter(vecM,vecMedianChange2,500,'k','_')

xlim([0 500]);
ylim([-100 110]);
yticks([-100 -50 0 50 100 110]);
yticklabels({'-100','-50','0','50','100','> 100'})
xlabel('Depth from surface (um)');
ylabel('Δ Rate (%)');
fixfig;
