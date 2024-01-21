%% plot the lick frequencies of mice, Hit and error separately, opto/no-opto in same plot
%script is written with output of analyseoptodata_lick_RT (b2r) as input

function [lickFrequencies, lickVecs4] = plotlickfrequencies(data,mousename_example)

%% prepare data summary (mouse x reaction x lickside x opto x timebin)
lickFrequencies = NaN(9,2,2,2,400);
mice = fieldnames(data);
numMice = length(mice);
stimLats = {'baseLats', 'oriLats', 'oopLats'};
stimReactions = {'baseReactions', 'oriReactions', 'oopReactions'};
stimRightVecs = {'baseRTrightVec', 'oriRTrightVec', 'oopRTrightVec'};
stimLeftVecs = {'baseRTleftVec', 'oriRTleftVec', 'oopRTleftVec'};
stimSides = {'baseSides', 'oriSides', 'oopSides'};
reactions = {'Hit', 'Error'};


for m = 1:numMice
    mn = mice{m};
    dat = data.(mn);
    gatherVecs = cell(2,2,2);
    
    for stim = 1:3
        if isfield(dat, stimLats{stim})
            trialsInclusion = ~(strcmp(dat.(stimReactions{stim}),'Miss')) & dat.(stimLats{stim}) <250;
            if sum(trialsInclusion) > 100
                optoVec = (dat.(stimLats{stim}) <250) + 1;
                for opto = 1:2 %no opto vs. opto
                    for reaction = 1:2 %hit vs. error
                        trialsSelect = find(strcmp(dat.(stimReactions{stim}),reactions(reaction)) & optoVec == opto);
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
    
    if strcmpi(mn, mousename_example)
        lickVecs4 = gatherVecs;
    end
    
    for opto = 1:2
        for reaction = 1:2
            for lickSide = 1:2
                if ~isempty(gatherVecs{opto,reaction,lickSide})
                    lickFrequencies(m,reaction,lickSide,opto,:) = mean(gatherVecs{opto,reaction,lickSide},1);
                else
                    lickFrequencies(m,reaction,lickSide,opto,:) = NaN;
                end
            end
        end
    end
end


end

