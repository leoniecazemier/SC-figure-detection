%% Lick data analysis script

loc = 'C:\Users\leoni\Documents\NIN\Data\Optomice';

%% get lick data and do statistics of RT
[b, goodses, PS] = analyseoptodata_lick_RT(loc, 0);

%% get summarized data and data from one mouse
[lickFrequencies, lickVecsMouse] = plotlickfrequencies(b,'Ibiza');

%% test max licking rates opto/no-opto (only hit, correct side)

%get mean frequencies per mouse, they are already smoothed
lickFrequencies_bylatency = getlickfrequencies_byoptolatency(b);
lickfreqs_hitcorrect = squeeze(lickFrequencies_bylatency(:,1,1,:,:));
lickfreqs_max = max(lickfreqs_hitcorrect,[],3);

latencies = [round([1:12]*(100/6)),250];

%remove mice without data
notnans = ~isnan(lickfreqs_max);
countnonans = sum(notnans,2);
lickfreqs_max = lickfreqs_max(countnonans>3,:);
%remove latencies without data
notnans = ~isnan(lickfreqs_max);
countnonans = sum(notnans,1);
lickfreqs_max = lickfreqs_max(:,countnonans > 3);
latencies = latencies(countnonans > 3);

%for anova
nmice = size(lickfreqs_max,1);
nlats = size(lickfreqs_max,2);

%use anovan 1 way rm anova - get stats
groupsM = 1:nmice;
groupsM = repmat(groupsM',1,nlats);
mousenums = reshape(groupsM,[],1);
groupsL = 1:nlats;
groupsL = repmat(groupsL',1,nmice);
latys = reshape(groupsL',[],1);
[PS.avnp,PS.avntbl,PS.avnstats] = anovan(reshape(lickfreqs_max,[],1),{latys,mousenums},'random',2);
PS.multcomp = multcompare(PS.avnstats);

% % bonferroni correction
% PS.OptoP = Pl* (size(lickfreqs_max,2)-1);

% Plot the results for max lick rate (of hit trials, licks on correct side)
figure(145)
clf
ph = plot(latencies,nanmean(lickfreqs_max,1),'o','MarkerSize',15,'MarkerFaceColor','k','MarkerEdgeColor', 'none','MarkerSize',15);
ebh = errorbar(latencies,nanmean(lickfreqs_max,1), nansem_large(lickfreqs_max,1),'.','color','k');
set(gca,'xtick',[0:50:250]);
set(gca, 'xticklabels', {'0', '50', '100', '150', '200','no'})
xlim([0 270])
ylim([5, 9])
xlabel('Laser onset (ms)')
ylabel('Max lick rate (licks/s)')
%

%% plot summarized data

px = -0.495:0.01:3.5;
ylims = [-8, 8; -8, 8];
optolabels = {'Laser OFF','laser ON'};
colors = {[0 0 0],[0.6 .6 .6];[0 0 1],[0 0.5 1]};
lickdirectionhelp = [1, -1];
reactions = {'Hit','Error','Miss'};

figure(141);
clf
plotcount = 0;
for reaction = 1:2 %hit/error
    for opto = 1:2 %off / on
        plotcount = plotcount+1;
        subplot(2,2,plotcount)
        cla;
        hold on;
        title([reactions{reaction}, ', ', optolabels{opto}])
        xlabel('Licks/s')
        xlabel('Time from stim onset (s)')
        yline(0);
        for lickside = 1:2 %correct/incorrect
            %correct licks, no opto
            [handlefill,handleline] = errorfill(px,smooth(nanmean(squeeze(lickFrequencies(:,reaction,lickside,opto,:)),1),10)*lickdirectionhelp(lickside), ...
                smooth(nansem_large(squeeze(lickFrequencies(:,reaction,lickside,opto,:)),1),10),colors{opto,lickside});
            alpha(handlefill,0.4)
        end
        ylim([-8, 8])
        xlim([-0.5, 2.5])
        ylabel('Licks/s')
    end
end


%% plot 1 trial

reactions = {'Hit','Error'};
optolabels = {'Laser OFF','laser ON'};
px = -0.495:0.01:3.5;
ylims = [-120, 120];
trial_no = 90;  %Aygo: 23
colors = {'-k','-b'};
lickdirectionhelp = [1, -1];

figure(215);
clf;
plotcount = 0;
for reaction = 1:2
    for opto = 1:2 % off/on
        plotcount = plotcount+1;
        subplot(2,2,plotcount)
        cla;
        hold on;
        title([reactions{reaction}, ', ', optolabels{opto}])
        xlabel('Licks/s')
        xlabel('Time from stim onset (s)')
        yline(0);
        
        for lickside = 1:2 %correct/incorrect
            
            licks{reaction,opto,lickside} = find(squeeze(lickVecsMouse{opto,reaction,lickside}(trial_no,:)) >50);
            if ~isempty(licks{reaction,opto,lickside})
                ax{reaction,opto,lickside} = plot([px(licks{reaction,opto,lickside}); px(licks{reaction,opto,lickside})], repmat([0 100*lickdirectionhelp(lickside)],1,size(licks{reaction,opto,lickside},1)), colors{opto}, 'linewidth', 2);
            end
        end
        ylim(ylims)
        xlim([-0.5, 2.5])
        ylabel('Licks')
        yticks([-100,100])
        yticklabels({'Incorrect','Correct'})
    end
end




