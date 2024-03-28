%% plot simulation results - Sec III methodology and sec IV Simulations
clear all; close all; clc


%% read files with simulated impulse responses
fs = 48000 * 1; % sampling rate - no oversampling by now

switch 2
    case 1
        for it = 1:2
            filename = 'IR_Synthetic_%d.wav'; % read IRs
            rir(:, it) = audioread(sprintf(filename,it));
        end
    case 2
        load('IR_ISM_small_change.mat');
        rir = h.ism;
        rir = rir(400:end, :);
end


len = length(rir);
noise = 0.00001*randn(len, 2);

minDB = 30;
%% bandpass filtering and correlation estimation
bandCenters = (1:20)*1000; % Hz - take only 1-19 kHz
win = rectwin(2^10);
win = win / sum(win);

for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt); % bandwidth 1000 Hz
    rir_band = bandpass(rir,band_freq,fs);
    noise_band = bandpass(noise,band_freq,fs);
    
    % compute correlation
    cor = slidingCorrelation(rir_band);
    [cor_noise, energy] = slidingCorrelation(rir_band + noise_band);

    meas_cor(:,bandIt) = cor(:,2);
    meas_cor_noise(:,bandIt) = cor_noise(:,2);

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:,bandIt) = energyDB > min(energyDB) + minDB ;
    mask(2*fs:end,bandIt) = 0;
   

    for it = 1:2 
       noise_energy_band(:, bandIt, it) =  conv(noise_band(:,it).^2,win,'valid');
       sig_energy_band(:, bandIt, it) =  conv(rir_band(:,it).^2,win,'valid');
    end
end
time_cor = (1:size(cor,1))/fs; % seconds

%% compute SNR-based correlation
med_en_sig = squeeze(median(sig_energy_band, 3));
med_en_noise = squeeze(median(noise_energy_band, 3));
snr_corr = med_en_sig./(med_en_sig + med_en_noise);
%% towards analytical model
[pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs);
div_ = 10^-5;
pred_time = pred_toastd / (sqrt(2)*div_); % double the divergence because of the difference of two Gaussians

%% check calibration towards baseline
% values at time t (1 second)
t = .2;
calibValuesSimu = interp1(time_cor, meas_cor, t)
calibValuesModel = interp1(pred_time, pred_cor, t)
%%
% load('div_sim.mat')
[divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask);
%% get the correlation coefficient between the calculated correlation curves and the fitted curves
t_cor = [0: length(cor)-1]./fs;
corL = length(meas_cor)./fs;
[~, loc] = min(abs(pred_time - corL));

new_cor = zeros(length(meas_cor), bandIt);
cor_PCC = zeros(2,2,bandIt);
for it = 1: bandIt
    new_cor(:, it) = interp1(pred_time(1:loc), pred_cor(1:loc, it), t_cor);
    cor_PCC (:, :, it) = corrcoef([meas_cor(:, it), new_cor(:, it)]);
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
%% plot figure 5 from the paper - measured vs single-pulse correlation of the simulated RIRs

% pred_time = 1000*pred_t;

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

figure(5); clf; hold on;
for i = 5
% plot(t_cor, meas_cor(:,i), 'color', cMap2(:, i), 'LineWidth', 1)
plot(t_cor, meas_cor_noise(:,i), '-', 'color', cMap2(:, i), 'LineWidth', 1)
plot(t_cor, meas_cor(:,i), 'r-.', 'LineWidth', 1)
plot(t_cor, snr_corr(1:length(t_cor),i), ':', 'color', cMap2(:, i), 'LineWidth', 1)
plot(pred_time,pred_cor(:,i),'--', 'color', cMap2(:, i), 'LineWidth', 1.5, 'HandleVisibility','off');
% plot(t_cor,new_cor(:,i).*snr_corr(:, i),'k--',  'LineWidth', 1.5, 'HandleVisibility','off');
end
xlabel('Time'); ylabel('Correlation');
legend('Simulated','Analytical');
xlim([0, 1.81])
ylim([-0.15, 1.01])

box on
set(gca, 'ytick', -0.25:0.25:1, 'FontSize', 12)
grid on

xlabel('Time (s)', 'Interpreter','latex')
ylabel('Correlation', 'Interpreter','latex')
% legend('$\rho_{h_i, h_j} $', '$\rho_{1\mathrm{r}} $', 'Interpreter', 'latex', 'FontSize',12)
% legend('5 kHz', '10 kHz', '19 kHz', 'Interpreter', 'latex', 'FontSize',12, 'location', 'southwest')
%% plot figure 6 from the paper - correlation at 1s, divergence vs reference value
figure(1); clf
colororder({'k', 'k'})
subplot(2,1,1); hold on
plot(bandCenters./1000, calibValuesSimu, 'k', 'LineWidth',1.5)
plot(bandCenters./1000, calibValuesModel, 'r-.', 'LineWidth',1.5)

xlim([1 20])
ylim([-0.01, 1.01])

box on
set(gca, 'ytick', -0.25:0.25:1, 'FontSize', 12)
grid on

xlabel('Center Frequecy (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation', 'Interpreter','latex')

legend('$\rho_{h_i, h_j} $', '$\rho_{1\mathrm{r}} $', 'Interpreter', 'latex', 'FontSize',12, 'location', 'southwest', 'numcolumns', 2)
subplot(2,1,2); hold on
plot(bandCenters./1000, divergence,'-k', 'LineWidth',1.5); hold on
plot(bandCenters./1000, ones(19,1)*10^-5,'r--',    'LineWidth',1.5)
xlim([1 19])
ylim([0.85 1.07]*10^-5)
labs = [0.85:0.05:1.05].*10^(-5);
set(gca, 'ytick', labs, 'YTickLabel' ,labs./10^(-5) )
text(0.85, 1.085*10^-5, '$\times 10^{-5}$', 'interpreter', 'latex', 'FontSize',12)
legend( '$\widetilde{\vartheta}$','Target', 'Interpreter', 'latex', 'FontSize',12, 'location', 'southwest', 'numcolumns', 2)
box on
grid on

xlabel('Center Frequecy (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('Divergence', 'Interpreter','latex')
set(gca,  'FontSize', 12)
% legend('$\widetilde{\vartheta}$', 'Interpreter', 'latex', 'FontSize',12)
%% figure 3 from the paper - modeled single reflection correlation over time and frequency

figure(6); clf; hold on;
for i = 1: bandIt
% plot(t_cor, meas_cor(:,i), 'color', cMap2(:, i), 'LineWidth', 1)
plot(pred_time,pred_cor(:,i),'-', 'color', cMap2(:, i), 'LineWidth', 1.5, 'HandleVisibility','off');
end

xlim([0, 2])
ylim([-0.02, 1.02])

box on
set(gca, 'ytick', -0.25:0.25:1, 'FontSize', 12)
grid on

xlabel('Time (s)', 'Interpreter','latex')
ylabel('Modeled correlation $\rho_{\textrm{1r}}$', 'Interpreter','latex')
% legend('$\rho_{h_i, h_j} $', '$\rho_{1\mathrm{r}} $', 'Interpreter', 'latex', 'FontSize',12)
clb = colorbar;
colormap(cMap2.')
clb.TickLabelInterpreter = 'latex';
clb.Ticks = .075:.1:1;
clb.TickLabels = 2:2:20;
clb.Label.String = 'Center frequency (kHz)';
clb.Label.Interpreter = 'Latex';
set(clb, 'Direction', 'reverse');
clb.FontSize = 12;
clb.Label.FontSize = 12;


%%
figure(7); clf; hold on
plot(bandCenters, squeeze(cor_PCC(1,2,:)), '-.')


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

pred_toastd = linspace(0,10^-3,10^5).'; % seconds; that's enough time for 1kHz to decay to 0
sigma = pred_toastd.' * fs; % toastd in samples

g = discreteGaussian(lags,sigma.^2);
% pred_cor = sum(c .* g, 1); % predicted correlation
pred_cor = g' * c; % predicted correlation

% makes a beautiful fit - Is twice the Gaussian equal to the power two?
% pred_cor = pred_cor.^0.5;

end


%% sliding correlation
function [xCor, energy] = slidingCorrelation(ir)
referenceIR = 1;
num_rirs = size(ir,2);
winLen = 2^12;
% win = ones(frameLen,1) / (frameLen);
win = rectwin(winLen);%hann(winLen); % much smoother than rectangular window
% win = ones(winLen,1);
win = win / sum(win);

for it = 1:num_rirs % for each channel
    xCov(:,it) = conv(ir(:,it) .* ir(:,referenceIR),win,'valid');
    energy(:,it) = conv(ir(:,it).^2,win,'valid');
    diffEnergy(:,it) = conv( (ir(:,it) - ir(:,referenceIR)).^2,win,'valid');
end
xCor = xCov ./ sqrt(energy) ./ sqrt(energy(:,referenceIR));

end

%% synthesize fractional delays
function rir = synthRIR(echo_times, echo_amps, len)

numberOfEchoes = size(echo_times,1);
RIR = zeros(len,1);

impulse = zeros(len,1);
impulse(2) = 1;

Impulse = fft(impulse);
A = 1i*angle(Impulse);

for it = 1:numberOfEchoes
    RIR = RIR + echo_amps(it) .* exp(A.* (echo_times(it)-1));
end

rir = real(ifft(RIR));

end

%% find divergence
function [divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask)
[numT, numIR, numBands] = size(meas_cor);
for itIR = 1:numIR
    for itBands = 1:numBands
        cor = meas_cor(:,itIR,itBands);
        
        m = mask(:,itBands);

        toastd_interp = interp1(pred_cor(:,itBands),pred_toastd,cor,'linear','extrap');

        b = robustfit(time_cor(m),toastd_interp(m),[],[],0);

        divergence(itIR,itBands) = b / sqrt(2); %median(localDivergence(mask(:,itBands)));

    end
end
end