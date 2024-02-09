%% Get data

load('N:\Ephys\Analysis\Singleunits\unitdatareviewed.mat')

%% compare depths of mp neurons with visual neurons

mpcs = find([unit.MP]);
vscs = find([unit.use] &([unit.useOri]|[unit.useOOP]) & ismember([unit.RFPOS],1:4) & [unit.actCell] & ~[unit.MP]);
[Val, Pos] = intersect(mpcs, vscs); %check for intersections

depths = [unit.depthR];

if isempty(Val)
    %get depths
    mpdepths = depths(mpcs);
    vsdepths = depths(vscs);
    
    %test for normality
    [Hmp, Pm] = lillietest(mpdepths);
    [Hvs, Pv] = lillietest(vsdepths);
    [Hvar,Pvar, CIvar, STvar] = vartest2(mpdepths, vsdepths)
    
    if Hmp || Hvs || Hvar
       [ Pd, Hd,Statsd] = ranksum(mpdepths,vsdepths);
    else
        [Hd, Pd] = ttest2(mpdepths, vsdepths);
    end
    
end


figure;
hold on;
VS = Violin(vsdepths,1,'ViolinColor',[0 0 0], 'ShowMean', true);
MP = Violin(mpdepths,2,'ViolinColor',[0 0 0], 'ShowMean', true);
xlim([0.5 2.5]);
ylim([0  750]);
set(gca, 'YDir','reverse')
set(gca, 'xtick', 1:2, 'linewidth', 1.2, 'fontsize', 14, 'Xticklabel', {'Visual','Multisensory'});
ylabel('Depth relative to SC surface (um)');

    
