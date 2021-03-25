%% Plot noisy histogram given noise scale
% Input:
%   noise_scale: scale of noise, add Lap(2 * noise_scale)
%   c: original Markov chain
%   plotoffset: offsite for plot position
%   dotcolor, marker: color/marker of the scatter plot
%   rep: number of repeatation 
% Return:
%   errl1: averaged l1 error
%   plt: handle of the plot 

function [errl1, plt] = plotHistNoise(noise_scale, c, plotoffset, dotcoloridx, marker, rep)

if nargin == 2
    plotoffset = 0;
    dotcoloridx = 6;
    marker = 'o';
    rep = 20;
end

color = get(groot,'DefaultAxesColorOrder');
if numel(dotcoloridx) == 1
    dotcolor = color(dotcoloridx,:);
else
    dotcolor = dotcoloridx;
end

rangeh = min(c):max(c);
h = hist(c, rangeh);
T = sum(h);
noise = gamrnd(1, 2*noise_scale, rep, length(h)) .* ((rand(rep, length(h)) > 0.5)*2-1);
h_noisy = repmat(h, rep, 1) + noise;

errl1 = mean(sum(abs(noise), 2))/T;

hold all
for j = 1:rep
    plt = scatter(rangeh+plotoffset, h_noisy(j,:)/T,50,'filled','Marker',marker, 'MarkerFaceAlpha',2/8,'MarkerFaceColor',dotcolor,'MarkerEdgeColor',dotcolor);
end


