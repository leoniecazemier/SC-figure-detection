%% plot the lick frequencies of mice, Hit and error separately, opto/no-opto in same plot
%script is written with output of analyseoptodata_lick_RT (b2r) as input

function [lickFrequencies_bylatency] = getlickfrequencies_byoptolatency(data)

%% prepare data summary (mouse x reaction x lickside x opto x timebin)
mice = fieldnames(data);
numMice = length(mice);
stimLats = {'baseLats', 'oriLats', 'oopLats'};
stimReactions = {'baseReactions', 'oriReactions', 'oopReactions'};
stimRightVecs = {'baseRTrightVec', 'oriRTrightVec', 'oopRTrightVec'};
stimLeftVecs = {'baseRTleftVec', 'oriRTleftVec', 'oopRTleftVec'};
stimSides = {'baseSides', 'oriSides', 'oopSides'};
reactions = {'Hit', 'Error'};
taskLats = [17, 33, 50, 67, 83, 100, 117, 133, 150, 167, 183, 200, 250];

lickFrequencies_bylatency = NaN(numMice,2,2,length(taskLats),400);


for m = 1:numMice
    mn = mice{m};
    dat = data.(mn);
    gatherVecs = cell(length(taskLats),2,2);
    
    for stim = 1:3
        if isfield(dat, stimLats{stim})
            trialsInclusion = ~(strcmp(dat.(stimReactions{stim}),'Miss')) & dat.(stimLats{stim}) <250;
            if sum(trialsInclusion) > 100
                lats_rounded = round(dat.(stimLats{stim}));
                for opto = 1:length(taskLats) %no opto vs. opto
                    optoVec = lats_rounded == taskLats(opto);
                    for reaction = 1:2 %hit vs. error
                        trialsSelect = find(strcmp(dat.(stimReactions{stim}),reactions(reaction)) & optoVec);
                        for trial = 1:length(trialsSelect)
                            tnum = trialsSelect(trial);
                            if strcmp(dat.(stimSides{stim})(tnum),'right')
                                RTvecCor = cell2mat(dat.(stimRightVecs{stim})(tnum));
                                RTvecIncor = cell2mat(dat.(stimLeftVecs{stim})(tnum));
                            else
                                RTvecCor = cell2mat(dat.(stimLeftVecs{stim})(tnum));
                                RTvecIncor = cell2mat(dat.(stimRightVecs{stim})(tnum));
                            end
                            RTidxCor = round((RTvecCor + 500)/10);
                            RTidxCor = RTidxCor(RTidxCor > 0 & RTidxCor<= 400);
                            saveVecCor = zeros(1,400);
                            saveVecCor(RTidxCor) = 100;
                            gatherVecs{opto,reaction,1} = [gatherVecs{opto,reaction,1}; saveVecCor]; %correct lick
                            RTidxIncor = round((RTvecIncor + 500)/10);
                            RTidxIncor = RTidxIncor(RTidxIncor > 0 & RTidxIncor<= 400);
                            saveVecIncor = zeros(1,400);
                            saveVecIncor(RTidxIncor) = 100;
                            gatherVecs{opto,reaction,2} = [gatherVecs{opto,reaction,2}; saveVecIncor]; %correct lick
                        end
                    end
                end
            end
        else
            continue
        end
    end
    
    for opto = 1:length(taskLats)
        for reaction = 1:2
            for lickSide = 1:2
                if ~isempty(gatherVecs{opto,reaction,lickSide})
                    temp_frequencies = mean(gatherVecs{opto,reaction,lickSide},1);
                    smoothed_frequencies = smooth(temp_frequencies,10);
                    lickFrequencies_bylatency(m,reaction,lickSide,opto,:) = smoothed_frequencies;
                else
                    lickFrequencies_bylatency(m,reaction,lickSide,opto,:) = NaN;
                end
            end
        end
    end
end


end

