%% plot the single reflection correlation against toa std and frequency
clear all; close all; clc

%%
fs = 48000;
%% bandpass filtering and correlation estimation
bandCenters = (1:20)*1000; % Hz
divergence = [0.5 1 2].*10^-5

%% towards analytical model
[pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs);

for i = 1:3
pred_time(:, i) = pred_toastd / (sqrt(2)*divergence(i)); % double the divergence because of the difference of two Gaussians


calibValuesModel(:, i) = interp1(pred_time(:, i), pred_cor, 1);
end
%% colors
numPlots = 20;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;

%% plots
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
f1 = figure(1); clf; hold on
x0=0;
y0=0;
width=80;
height=60;
% set(gcf,'position',f1.Position + [x0,y0,width,height])

% sf1 = subplot(2,1,1)
hold on
for i = 1:19
plot(pred_toastd*1000, pred_cor(:, i), 'color', cMap2(:, i), 'LineWidth',1)
end
xlim([0 0.1])
ylim([-0.01 1.01])
box on
grid on
set(gca, 'YTick',0:.2:1, 'FontSize',12)
xlabel('$\textrm{TOASTD}~\tau_\sigma$ (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('$\textrm{Single-reflection correlation}~\rho_{\textrm{1r}}$', 'Interpreter','latex', 'FontSize',12)
clb = colorbar;
colormap(cMap2.')
clb.TickLabelInterpreter = 'latex'
clb.Ticks = .075:.1:1
clb.TickLabels = 2:2:20
clb.Label.String = 'Center frequency (kHz)'
clb.Label.Interpreter = 'Latex'
set(clb, 'Direction', 'reverse')
clb.FontSize = 12
%%

% sf2 = subplot(2,1,2)
% sf2.Position(3) = sf1.Position(3)-0.008;
figure(2); clf
hold on
for i = 1:3
plot(bandCenters./1000, calibValuesModel(:, i), 'color', cMap2(:, 9*(i-1)+1), 'LineWidth',1.5)
end
xlim([1 19])
ylim([-0.01 1.01])
set(gca, 'YTick',0:.2:1, 'FontSize',12)
box on
grid on
xlabel('Center Frequecy (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('$\textrm{Single-reflection correlation}~\rho_{\textrm{1r}}$', 'Interpreter','latex', 'FontSize',12)
lgd = legend('0.005', '0.01','$0.02$', 'interpreter', 'latex', 'location', 'southwest', 'fontsize',12);
lgd.Title.String = '$\tau_\sigma$ (ms)'
% lgd.Position(1) = lgd.Position(1)+0.15

%% half life for different divergences
divergence_ = [0.5:0.1:2].*10^-5

%% towards analytical model
[pred_cor_,pred_toastd_] = singlePulseCorrelation(bandCenters,fs);
const = 1/exp(1);
[~, loc] = min(abs(pred_cor_-const),[],1);

for i = 1:length(divergence_)    
    pred_time_(:, i) = pred_toastd_ / (sqrt(2)*divergence_(i)); % double the divergence because of the difference of two Gaussian
    half_time(:, i) = pred_time_(loc, i);
% calibValuesModel_(:, i) = interp1(pred_time_(:, i), pred_cor_, 1);
end


%% carbon half life of correlation
figure(3); clf;
hold on
for i = 1:19
plot(divergence_, half_time(i,:), 'color', cMap2(:, i), 'LineWidth',1)
end 

% xlim([0 0.1])
ylim([0.35 35])
box on
grid off
grid on
set(gca, 'YTick',[0.5 1 2 5 10 15 25], 'FontSize',12)
xlabel('Divergence $\vartheta$ (s)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation time (s)', 'Interpreter','latex', 'FontSize',12)
set(gca, 'yscale', 'log')
clb = colorbar;
colormap(cMap2.')
clb.TickLabelInterpreter = 'latex'
clb.Ticks = .075:.1:1
clb.TickLabels = 2:2:20
clb.Label.String = 'Center frequency (kHz)'
clb.Label.Interpreter = 'Latex'
set(clb, 'Direction', 'reverse')
clb.FontSize = 12
%% single pulse correlation
function [pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs)

fs = fs*8; % oversample filters, important for the highest bands

for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);

    pulse = zeros(fs,1); pulse(end/2) = 1;
    pulse_band = bandpass(pulse,band_freq,fs);

    max_lag = fs / 100; % this needs to be long enough for the filter
    [c(:,bandIt),lags] = xcorr(pulse_band,max_lag,'normalized'); % measure the autocorrelation of a bandpassed pulse
    lags = lags.'; % samples
end

pred_toastd = linspace(0,10^-3,10^4).'; % seconds; that's enough time for 1kHz to decay to 0
sigma = pred_toastd.' * fs; % toastd in samples

g = discreteGaussian(lags,sigma.^2);
% pred_cor = sum(c .* g, 1); % predicted correlation
pred_cor = g' * c; % predicted correlation

% makes a beautiful fit - Is twice the Gaussian equal to the power two?
% pred_cor = pred_cor.^0.5;

end
