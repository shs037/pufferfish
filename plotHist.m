%% Plot non-private histogram from original Markov chain

function plt = plotHist(c)
rangeh = min(c):max(c);
h = hist(c, rangeh);
T = sum(h);
plt = bar(rangeh, h/T, 'facecolor', .9*ones(1,3),'EdgeColor',zeros(1,3));


