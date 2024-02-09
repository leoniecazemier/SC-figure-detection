%%% Plot eye movement
%% Get data
% load 'N:\Ephys\Analysis\behTanks\eyeDat_Z.mat'
% load 'N:\Ephys\Analysis\Singleunits\unitdatareviewed_trimmed.mat'
load 'C:\Users\leoni\Documents\NIN\Data\eyeDat_Z.mat'
load 'C:\Users\leoni\Documents\NIN\Data\unitdatareviewed_trimmed.mat'

%% Settings
px = -0.1995:0.001:1.7495;
pxl =  -1.4995:0.001:0.9995;
sessions = unique({unit.ses});

% timebins
px = -0.1995:0.001:1.7495;
pxl = -1.4995:0.001:0.9995;

% Some info
trtypes = [13 14 15 16 17 18; 1 2 3 4 5 6; 7 8 9 10 11 12];  %figgrey grgrey figori grori figoop groop
fgTypeSplit = [1 2 3; 4 5 6];
sesList = {unit.ses};
seslist_unique = unique(sesList);
mice_names = unique({unit.mouse});
n_mice = length(mice_names);

%cell selection
rfSel = 1; %center
cSelection = [unit.use] & ismember([unit.RFPOS],rfSel) & [unit.actCell];
oriCellSel = [unit.useOri] & cSelection;
oopCellSel = [unit.useOOP] & cSelection;
conCellSel = ([unit.useOri]|[unit.useOOP]) & cSelection;
taskCellSel = {conCellSel oriCellSel oopCellSel};

%for plotting
colorsFG = {[0 210 75]./256,[0 102 51]./256,[235 32 39]./256,[121 10 10]./256,[68 150 200]./256,[31 67 132]./256 };
colorsP = {[190, 112, 250]/256,[95, 50, 120]/256};
rfNames = {'Inside','Outside', 'Edge', 'Unclear','All'};
fgNames = {'Figure', 'Ground'};
taskNames = {'Contrast','Orientation','Phase'};
eye_parameter_names = {'Eye X position','Eye Y position','Pupil dilation'};
heNames = {'Hit', 'Error'};
smfvis = 20;
smflick = 50;
plotIdxHelp = [2 3 1];

%some sessions it appears one paramater of eye data exists out of noise, probably cable
%malfunction. also first two mice had no eye data. Register which ones to use
eye_data_available = NaN(20,3);
eye_data_available(1:12,:)              = repmat([false, false, false],12,1);  %xpos, ypos, dilation
eye_data_available(13:14,:)             = repmat([false, true, true],2,1);  %xpos, ypos, dilation
eye_data_available([16, 18, 19, 20],:)  = repmat([true, false, true],4,1);  %xpos, ypos, dilation
eye_data_available([15,17],:)           = repmat([true, true, true],2,1);  %xpos, ypos, dilation


%% plot raw trials for one mouse

fieldNames = {'xpos_Z','ypos_Z','dil_Z','xposL_Z','yposL_Z','dilL_Z'};
axisNames = {'Eye X pos (Z)','Eye Y pos (Z)','Pupil dilation (Z)'};
session = 15;
sessionName = sessions{session};
smfvis = 20;
plot_units = find(strcmp({unit.ses}, sessionName) & [unit.use]) ;
clear traces
for trial = 59
    figure(trial);
    clf;
    for isubplot=1:3
        subplot(4,1,isubplot)
        plot(px,eyeDat(session).(fieldNames{isubplot})(trial,:));
        traces(isubplot,:) = eyeDat(session).(fieldNames{isubplot})(trial,:);
        xlabel('Time from stim onset (s)')
        ylabel(axisNames{isubplot});
    end
    subplot(4,1,4)
    clear resp
    for i_unit = 1:length(plot_units)
        c_unit = plot_units(i_unit);
        resp(i_unit,:) = unit(c_unit).nrDataVis(trial,:);
    end
    plot(px,nanmean(resp,1));
    traces(4,:) = nanmean(resp,1);
    xlabel('Time from stim onset (s)')
    ylabel('Normalized firing rate');
end



%% Plot eye position and pupil size for each task, f vs. g, across mice
eye_data = cell(n_mice,3,2); %mouse x xyd x fg
figure(21)
clf;


for mouse = 1:n_mice
    sessions_mouse = unique(sesList(strcmp({unit.mouse},mice_names(mouse))));
    mouse_task_data = cell(3,2);  %xyd x fg
    
    for session = 1:length(sessions_mouse)
        unit_options = find(strcmp(sesList,sessions_mouse(session)));
        unit_id = unit_options(1);
        has_eye_data = unit(unit_id).eyedata;
        if ~ has_eye_data
            continue
        end
        session_idx = find(strcmp(seslist_unique,sessions_mouse(session)));
        eye_data_ses = eyeDat(session_idx);
        trial_groups = unit(unit_id).trialGroups;
        %get traces and remove trials with artefacts (> 7 SD)
        for fg = 1:2
            trials_use = find(ismember(trial_groups, trtypes(:,[1:3] + (3*(fg-1)))));
            if eye_data_available(session_idx,1)
                temp = [mouse_task_data{1,fg}; eye_data_ses.xpos_Z(trials_use,:)];
                mouse_task_data{1,fg} = temp(max(abs(temp),[],2)<5,:);
            end
            if eye_data_available(session_idx,2)
                temp = [mouse_task_data{2,fg}; eye_data_ses.ypos_Z(trials_use,:)];
                mouse_task_data{2,fg} = temp(max(abs(temp),[],2)<5,:);
            end
            if eye_data_available(session_idx,3)
                temp = [mouse_task_data{3,fg}; eye_data_ses.dil_Z(trials_use,:)];
                mouse_task_data{3,fg} = temp(max(abs(temp),[],2)<5,:);
            end
        end
        
    end
    if all(cellfun(@isempty,mouse_task_data),[1,2])
        continue
    end
    for eye_parameter = 1:3
        for fg = 1:2
            if ~isempty(mouse_task_data{eye_parameter,fg})
                eye_data{mouse,eye_parameter,fg} =  nanmean(mouse_task_data{eye_parameter,fg},1);
            end
        end
    end
end

for eye_parameter = 1:3
    n_plot = eye_parameter;
    subplot(3,1,n_plot)
    hold on;
    clear mean_eye_data sem_eye_data
    %get mean and sem
    prep_fig = cat(1,eye_data{:,eye_parameter,1});  %concat across mice
    prep_gnd = cat(1,eye_data{:,eye_parameter,2});
    mean_eye_data(1,:) = smooth(nanmean(prep_fig,1),10); %fig
    mean_eye_data(2,:) = smooth(nanmean(prep_gnd,1),10); %gr
    sem_eye_data(1,:) = smooth(nansem_large(prep_fig,1),10); %fig
    sem_eye_data(2,:) = smooth(nansem_large(prep_gnd,1),10); %gr
    color_f = [190, 112, 250]/256;
    color_g = [95, 50, 120]/256;
    [handlefill(1), handleline(1)] = errorfill(px,mean_eye_data(1,:),sem_eye_data(1,:),color_f);
    [handlefill(2), handleline(2)] = errorfill(px,mean_eye_data(2,:),sem_eye_data(2,:),color_g);
    title(eye_parameter_names{eye_parameter});
    %         xlim([-0.05 0.25])
    xlim([-0.2 1.5])
    if eye_parameter == 1
        ylim([-0.6, 0.6])
        legend(handleline,{'Figure','Ground'},'location', 'northwest');
    elseif eye_parameter ==2
        ylim([-1.3, 1.3])
    else
        ylim([-0.5, 1.25])
    end
    
    xticks([0 0.5 1 1.5])
    
    
end

