function [behaviourTemp] = getoptoresults_3r(dat,useTr,mn,trtype, plotit)

LEG = {'Hit', 'Error','Miss'};
Green = [0 0.8 0];
Red = [0.8 0 0];
Grey = [0.5 0.5 0.5];

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
    P(l) = H(l)/(H(l)+E(l)+M(l));
    N(l) = H(l)+E(l)+M(l);
    MP(l) = M(l)/length(pickreact);
end

pickreact = react(isnan(lats));
H_NO = sum(strcmp(pickreact, 'Hit'));
E_NO = sum(strcmp(pickreact, 'Error'));
M_NO = sum(strcmp(pickreact, 'Miss'));
P_NO = H_NO/(H_NO+E_NO+M_NO);
MP_NO = M_NO/length(pickreact);

H = [H H_NO];
E = [E E_NO];
M = [M M_NO];
P = [P P_NO];
MP = [MP MP_NO];
mylats = [ulats 0.25];

behaviourTemp = {H; E; P; mylats; M};
toPl = [H; E; M];
mySum = sum(toPl,1);
toPl = toPl./mySum;

if plotit
    figure('Position', [65  200 1400   500]);
    subplot(1,2,1)
    mybar = bar(mylats,toPl','stacked');
    set(mybar,{'FaceColor'},{Green;Red;Grey});
    ylim([0 1]);
    title([trtype ' performances ' mn])
    ylabel('Fraction correct responses');
    for ll = 1:length(mylats)-1
        xval = mybar.XData;
        offs = mybar.XOffset;
        text(xval(ll)-0.005, 0.96, num2str(N(ll)))
    end
    xticks([mylats(1:end-1) 0.25]);
%     line([0 0.3], [0.5 0.5],'color', 'g', 'linestyle','-')
    line([0 0.3], [P_NO P_NO],'color', 'r', 'linestyle','-')
    labelhelp = arrayfun(@(x) num2str(x,'%3.0f'), mylats(1:end-1)*1000,'un',0);
    xticklabels([labelhelp 'No opto']) ;
    legend(LEG, 'location','southoutside')
    
end


end
