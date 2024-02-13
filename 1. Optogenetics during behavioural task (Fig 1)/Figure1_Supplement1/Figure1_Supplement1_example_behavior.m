%Figure1_Supplement1_example_behavior
%
% 2023, Alexander Heimel

%% Parameters and folders
rootpath = '.';
project = '.';

%% Constants
types = {'Gray','Even','Low contrast',...
    'Full contrast'};
n_types = length(types);

cm = parula(n_types+4);
cm = cm(4:end-1,:);
% contrast = 0x00d14b = 0,209,75
% orient = 0xea2027 = 234,32,39

cm(1,:) = [  0 209  75]/255; % Gray background
cm(2,:) = [150 180  63]/255; % Even background
cm(3,:) = [200 120  51]/255; % Low contrast
cm(4,:) = [234  32  39]/255; % Full contrast

%Source => Target = (BGColor + Source) =
%Target.R = ((1 - Source.A) * BGColor.R) + (Source.A * Source.R)
%Target.G = ((1 - Source.A) * BGColor.G) + (Source.A * Source.G)
%Target.B = ((1 - Source.A) * BGColor.B) + (Source.A * Source.B)
% 
% going to alpha = 0.8 (as preferred by Leonie)
alpha = 0.8;
cm = (1-alpha)*1 + alpha*cm;


laser_color = [41,170,225]/255; % 0x29aae1

% Opto mice in manuscript: Aygo, Elantra, Fiesta, Ibiza, Kangoo, Leaf, Mustang, Niro. 

dataset = 'Behaviour - opto mice';
mouse = 'Ibiza';
session_nr = 59;   % alternatives 80, 70, 67

data_path = fullfile(rootpath,project,'Data_collection',dataset,mouse);

%% Make single session figure

% Load session data
d = dir(fullfile(data_path,'*B1.mat'));
load(fullfile(data_path,d(session_nr).name));
disp(['Loaded ' d(session_nr).name])

% Extract data
correct_is_left = strcmp(LOG.Side,'left');
correct_is_right = strcmp(LOG.Side,'right');
hits = strcmp(LOG.Reaction,'Hit') & ~LOG.Gavepassive;
misses = strcmp(LOG.Reaction,'Miss') & ~LOG.Gavepassive;
errors = strcmp(LOG.Reaction,'Error') & ~LOG.Gavepassive;
reactions = ones(size(hits));
reactions(hits) = 3;
reactions(errors) = 2;

% CurrPerformance is sum(hits)|sum(hits|errors) over last 15 hits and errors

% Make performance figure
figure
hold on
plot(LOG.CurrPerformance ,'ok')
plot(hits,'.g' );
plot(errors,'.r' );
xlabel('Trial')
session_performance = sum(hits)/sum(hits|errors);
disp(['Session performance: ' num2str(round(session_performance*100)) ' %'])
plot(cumsum(hits)./cumsum(hits|errors),'b-')

%% Figure 1-figure supplement 1A. figure with ephys session behavior lick
figure;
n_trials = length(LOG.Trial);
hold on
if isfield(LOG,'optotrial')
    hb = barh(1:n_trials,2.6*LOG.optotrial,'BarWidth',1,'FaceColor',laser_color,'FaceAlpha',1,'EdgeColor','none');
    disp(['Number of optotrials: ' num2str(sum(LOG.optotrial))]);
    disp(['Performance on optotrials: ' num2str(round(sum(LOG.optotrial & hits) / sum(LOG.optotrial & (hits | errors))*100)) ' %']);
end
image('xdata',[0 1],'ydata',[1 n_trials],'cdata',reactions')
lick_colormap = [1 1 1;0.5 0.5 0.5;0.8 0.8 0.8];
colormap(lick_colormap)
for i = 1:n_trials
    r = LOG.RTleftVec{i};
    r(r>2000) = [];
    if ~isempty(r)
        hl = plot(r/1000,i*ones(size(r)),'r.');
    end
    r = LOG.RTrightVec{i};
    r(r>2000) = [];
    if ~isempty(r)
        hr = plot(r/1000,i*ones(size(r)),'g.');
    end
end % i


xlim([0 2.05])
ylabel('Trial #');
xlabel('Time from stim onset (s)')
hleg = legend([hl,hr,hb],{'Left lick','Right lick','Opto trial'},'location','northeastoutside','fontsize',12);
legend box off
axis square
legend_position = get(hleg,'position');
cb = colorbar;
set(cb,'ticks',[1.33 2 2.66]);
set(cb,'ticklabels',{'Miss','Error','Hit'})
set(cb,'FontSize',12)
p = get(cb,'position');
set(cb,'position',[legend_position(1)+0.05,p(2),0.05,legend_position(2)-p(2)])
ylim([0 n_trials])
%ylim([0 190])

%% Figure 1-figure supplement 1B. Make all session performance figure
data_path = fullfile(rootpath,project,'Data_collection',dataset,mouse);

d = dir(fullfile(data_path,[mouse '*B1.mat']));
dates = arrayfun( @(x) x.name(length(mouse)+2:length(mouse)+9),d,'UniformOutput',false);
dates = unique(dates);
n_sessions = length(dates);

session_performances = NaN(1,n_sessions);
training_types = {};
opto_sessions = false(1,n_sessions);
for i = 1:n_sessions
    d = dir(fullfile(data_path,[mouse '_' dates{i} '_*B*.mat']));
    hits = [];
    misses = [];
    errors = [];
    for j = 1:length(d)
        load(fullfile(d(j).folder,d(j).name));
        if ~isfield(LOG,'optotrial')
            LOG.optotrial = false(size(LOG.Reaction));
        end
        hits = [hits strcmp(LOG.Reaction,'Hit')& ~LOG.Gavepassive & ~LOG.optotrial]; %#ok<*AGROW>
        misses = [misses strcmp(LOG.Reaction,'Miss')& ~LOG.Gavepassive & ~LOG.optotrial];
        errors = [errors strcmp(LOG.Reaction,'Error')& ~LOG.Gavepassive & ~LOG.optotrial];
        if isfield(LOG,'TrainingType')
            training_types{i} = LOG.TrainingType; %#ok<*SAGROW>
            switch training_types{i}
                case {'Contrast','BaseTask','Texture'}
                    training_types{i} = 'Gray';
                    if any(LOG.BGContrast>0)
                        training_types{i} = 'Low contrast'; % not sure if this is correct
                    end
                case 'Luminance'
                    training_types{i} = 'Even';
                case 'Contrast Texture'
                    training_types{i} = 'Low contrast';
                case {'Full','Combi','Full Texture'} 
                    training_types{i} = 'Full contrast';
            end
        else
            training_types{i} = 'Gray'; % not sure if this is correct
        end
        if isfield(LOG,'optotrial') && any(LOG.optotrial)
            opto_sessions(i) = true;
        end
    end
    if sum(hits|errors)<2
        session_performances(i) = NaN;
    else
        session_performances(i) = sum(hits)/sum(hits|errors);
    end
end

session_types = zeros(1,n_sessions);
for i=1:length(types)
    session_types(strcmp(training_types,types{i}))=i;
end

figure
hold on
bar((1:n_sessions),110*opto_sessions,'facecolor',laser_color,'facealpha',1,'barwidth',1,'edgecolor','none');
image('xdata',1:n_sessions,'ydata',[0 50],'cdata',session_types);
colormap(cm)
stairs((1:n_sessions+1)-0.5,[session_performances session_performances(end)]*100,'-k')
ylim([0 110])
xlabel('Session');
ylabel('Performance (%)')
plot(xlim,[50 50],'-','color',[0.7 0.7 0.7])
cb = colorbar();
set(cb,'ticks',0.6+(1:n_types)*0.75);
set(cb,'ticklabels',types)
set(cb,'FontSize',12)
set(cb,'Direction','reverse')
p = get(gca,'position');
pcb = get(cb,'position');
set(cb,'position',[pcb(1) 0.5 pcb(3) pcb(4)+pcb(2)-0.5])
set(cb.Label,'String','Background')
set(gca,'position',p)
bar((1:n_sessions),110*opto_sessions,'facecolor',laser_color,'facealpha',0.2,'barwidth',1,'edgecolor','none');
xlim([1 length(session_performances)+0.5])
plot(session_nr,session_performances(session_nr)*100+5,'v','color',[0 0 1],'markerfacecolor',[0 0 1])