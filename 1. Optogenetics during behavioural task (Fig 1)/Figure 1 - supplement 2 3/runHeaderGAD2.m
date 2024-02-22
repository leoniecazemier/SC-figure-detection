
%define paths
addpath(fullfile(fileparts(mfilename('fullpath')),'dependencies'));

strSavePath = fullfile(fileparts(mfilename('fullpath')),'preprocessed_data');

if ~exist(strDataPath,'dir')
    warning('Data path does not exist!');
    sExp = [];
    return
end
sFiles = dir(fullpath(strDataPath,'*_AP.mat'));
if ~exist('sExp','var') || isempty(sExp)
    sExp = [];
    for intFile=1:numel(sFiles)
        fprintf('Loading %d/%d: %s [%s]\n',intFile,numel(sFiles),sFiles(intFile).name,getTime);
        sLoad = load(fullpath(sFiles(intFile).folder,sFiles(intFile).name));
        sLoad.sAP.Name = sFiles(intFile).name;
        if isempty(sExp)
            sExp = sLoad.sAP;
        else
            sExp(end+1) = sLoad.sAP; %#ok<SAGROW>
        end
    end
end

cellUseAreas = {'Superior colliculus zonal layer', 'Superior colliculus superficial gray layer', 'Superior colliculus optic layer'};
cellSubjectGroups = {'GAD2Cre','BL6'};
