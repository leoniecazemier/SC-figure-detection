function clusters = ez_clusterstat_time(cond1,cond2,mice, sessions, units, figon,t)

%[sigmap,threshmap,pixmap] = ez_clusterstat(cond1,cond2,freq,t,figon)
%Carries out cluster-based staistical test for paired data as described by
%Maris and Oostenveld (2007). Designed to work with baseline-corrected
%time-based data which has only positive clusters (i.e. the test is
%one-tailed, H1: cond1 is larger than cond2)

%The scrambling occurs by switching condition labels for each electrode.
%Electrode pairings are maintained (the data is assumed to be paired).

%Matt Self 2022

%Inputs:
%cond1 and cond2: 2d matrices with dimensions [time, electrode]
%Any missing data should be converted to NaNs before running this script
%These are the two conditions you want to compare e.g. figure and ground

%figon = draws figures

%Outputs:
%clusters = structure containing logical map of significant time points for
%each significant cluster.


%APPROACH
%1: MAke a permutation surrogate of the data by randomly shuffling conditions for
%each electrode to mak
%2: Calculate the t- and p-statistic between cond1' and
%   cond2' using a paired t-test for every sample.
%3: Cluster the fmap after thresholding with the pmap
%4: Sum-up the absolute t-vals from each cluster to get a distribution
%   of summed cluster t-values.  We take the maximum cluster t-value
%   as the bootstrap statistic as this controls for multiple comparisons.
%5: Do this 1000 times to work out the critical cluster maximum t-value for the 5% familywise alpha
%   level
%6: Apply the same procedure to the real data, clusters falling above the
%   critical value are significant.
%7: Mark these clusters in a graph

if nargin<6
    figon = 0;
    %Only used for figures
    t = [];
end

%% data dimensions
nsamps = size(cond1,1);
nelectrodes = size(cond1,2);

%% MAke a joint distribution contiang the data from both conditions
J = cat(3,cond1,cond2);
    
%% Number of permutations, 1000 is reasonable
reps = 1000;

%% perumtation stats
fmax = zeros(reps,1);

conditions = [repmat(1,1,size(cond1,2)),repmat(2,1,size(cond1,2))];
mice_all = repmat(mice,1,2);
sessions_all = repmat(sessions,1,2);
units_all = repmat(units,1, 2);


for s = 1:reps
    
    if rem(s-1,100) == 0
        fprintf('Bootstrapping: working on repetitions %i until %i \n',s, s+99);
    end

    %Randomly permute the conditions
    %E.g. either swap or do not swap the conditions for each electrode
    c = randi(2,nelectrodes,1);
    
    %Draw a new set of conditions from the joint distribution
    c1 = zeros(size(cond1));
    c2 = zeros(size(cond2));
    for n = 1:nelectrodes
        c1(:,n) = squeeze(J(:,n,c(n)));
        c2(:,n) = squeeze(J(:,n,3-c(n)));
    end
    
    %PErform the t-test
    c_all = [c1 c2];
    pmap = nan(1,size(c1,1));
    fmap = nan(1,size(c1,1));
    for tbin = 1:size(c1,1)
        [pmap(tbin),fmap(tbin)] = pairedNestedDataPValue(c_all(tbin,:)',conditions', mice_all',sessions_all', units_all');
    end
%     tmap = stats.tstat;
    
    %Threshold with the pmap into positive clusters
    fmap_pos = zeros(size(fmap));
    fmap_pos(fmap>0&pmap<0.05) = fmap(fmap>0&pmap<0.05);

    %Perform clustering
    %get labeled map via bwconncomp
    blobinfo = bwconncomp(fmap_pos);
    nblobs = blobinfo.NumObjects;
    if nblobs>0
        clustsum   = zeros(1,nblobs);
        for i=1:nblobs
            clustsum(i) = sum(fmap(blobinfo.PixelIdxList{i}));
        end
        fmax(s,1) = max(abs(clustsum));
    end
end


%% Calculate the maximum cluster statistic
J = sort(fmax);
%95th percentile for one-tailed test
percentile = round(reps.*0.95);
%This is the critical vlaue of maximum t
cluscrit = J(percentile);

%% Now cluster the real data and check to see whether each cluster was significant
%PErform the t-test
cond_all = [cond1 cond2];
pmap = nan(1,size(cond1,1));
fmap = nan(1,size(cond1,1));
for tbin = 1:size(cond1,1)
    [pmap(tbin),fmap(tbin)] = pairedNestedDataPValue(cond_all(tbin,:)',conditions', mice_all',sessions_all', units_all');
end
%% Cluster level correction
%Threshold with the pmap
fmap_pos = zeros(size(fmap));
fmap_pos(fmap>0&pmap<0.05) = fmap(fmap>0&pmap<0.05);

%Intiialise outputs
clusters = [];
cc = 0;

% get labeled map
blobinfo = bwconncomp(fmap_pos);
nblobs = blobinfo.NumObjects;
if nblobs>0
    clustsum   = zeros(1,nblobs);
    for i=1:nblobs
        clustsum(i)   = sum(fmap(blobinfo.PixelIdxList{i}));
    end
    
    %Clusters larger than threshold
    clustix = find(abs(clustsum)>cluscrit);
    if ~isempty(clustsum)
        for j = clustix
            buf = zeros(nsamps,1);
            buf(blobinfo.PixelIdxList{j}) = 1;
            cc = cc+1;
            clusters(cc).map = buf;
        end
    end
end


if figon
    
    figure
    plot(t,log10(pmap)),hold on
    title('Uncorrected stats'),ylabel('log(P)')
    
    %Mark significant clusters on map
    figure
    conddiff = nanmean(cond1-cond2,2);
    plot(t,conddiff,'k','LineWidth',2),hold on
    Y = get(gca,'YLim');
    for j = 1:length(clusters)
        st = find(clusters.map,1,'first');
        ed = find(clusters.map,1,'last');
        fill([t(st) t(st) t(ed) t(ed)],[Y(1) Y(2) Y(2) Y(1)],'g')
    end
    hold on,plot(t,conddiff,'k','LineWidth',2)
    title('Significant Clusters')
    
end




%% Colormap generation
function semap = makesemap(minval,maxval)

%MAke a suppression/enhancement colormap
cspace = linspace(minval,maxval,64);

%Find zero-point
zp = find(cspace>=0,1,'first');
semap = zeros(64,3)+0.5;
semap(zp:end,1) = linspace(0.5,1,64-zp+1)';
semap(zp:end,2) = linspace(0.5,0,64-zp+1)';
semap(zp:end,3) = linspace(0.5,0,64-zp+1)';
%All pojts below scale from grey to blue
semap(1:zp,1) = linspace(0,0.5,zp)';
semap(1:zp,2) = linspace(0,0.5,zp)';
semap(1:zp,3) = linspace(1,0.5,zp)';

return



