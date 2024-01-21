function [behaviourTemp] = getoptoresults_2r(dat,useTr,mn,trtype, plotit)


lats = dat.optolats(useTr);
react = dat.reactions(useTr);
ulats = unique(lats(~isnan(lats)));
clear H E P MP
for l = 1:length(ulats)
    lat = ulats(l);
    pickreact = react(lats == lat);
    H(l) = sum(strcmp(pickreact, 'Hit'));
    E(l) = sum(strcmp(pickreact, 'Error'));
    M(l) = sum(strcmp(pickreact, 'Miss'));
    P(l) = H(l)/(H(l)+E(l));
    N(l) = H(l)+E(l);
    MP(l) = M(l)/length(pickreact);
end

pickreact = react(isnan(lats));
H_NO = sum(strcmp(pickreact, 'Hit'));
E_NO = sum(strcmp(pickreact, 'Error'));
M_NO = sum(strcmp(pickreact, 'Miss'));
P_NO = H_NO/(H_NO+E_NO);
MP_NO = M_NO/length(pickreact);

P = [P P_NO];
H = [H H_NO];
E = [E E_NO];
MP = [MP MP_NO];
mylats = [ulats 0.25];

behaviourTemp = {H; E; P; mylats};

if plotit
    figure('Position', [65  200 1400   500]);
    subplot(1,3,1)
    mybar = bar(mylats,P);
    ylim([0 1]);
    title([trtype ' performances ' mn])
    ylabel('Fraction correct responses');
    for ll = 1:length(mylats)-1
        xval = mybar.XData;
        offs = mybar.XOffset;
        text(xval(ll)-0.005, 0.96, num2str(N(ll)))
    end
    xticks([mylats(1:end-1) 0.25]);
    line([0 0.3], [0.5 0.5],'color', 'g', 'linestyle','-')
    line([0 0.3], [P_NO P_NO],'color', 'r', 'linestyle','-')
    labelhelp = arrayfun(@(x) num2str(x,'%3.0f'), mylats(1:end-1)*1000,'un',0);
    xticklabels([labelhelp 'No opto']) ;
    
    subplot(1,3,2)
    mybar = bar(mylats,MP);
    ylim([0 1]);
    title([trtype ' perc misses ' mn])
    xticks([mylats(1:end-1) 0.25]);
    labelhelp = arrayfun(@(x) num2str(x,'%3.0f'), mylats(1:end-1)*1000,'un',0);
    xticklabels([labelhelp 'No opto']) ;
    line([0 0.3], [MP_NO MP_NO],'color', 'r', 'linestyle','-')
    ylabel('Fraction missed trials');
end

end
