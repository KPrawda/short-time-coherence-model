% Sebastian J. Schlecht, Friday, 22. December 2023
% test SNR correlation
clear; clc; close all;

% parameters
fs = 48000;
len = fs;
winLen = 400;
noiseLevel = 0.1;
rirLevel = linspace(1,0.2,len).';
sinv = @(snr) clip(snr ./ (1 + snr),[eps 1]); % snr to correlation conversion

% synthesize signal
noise = randn(len,2) .* noiseLevel;
rir = randn(len,2) .* rirLevel;

% sparsify signal
% rir(randperm(len,len/2),:) = eps;

rir_correlation = 0.95;
target_correlation = [1, rir_correlation; rir_correlation, 1];
mixing_matrix = chol(target_correlation);
rir = rir * mixing_matrix;
xcorr(rir(:,1),rir(:,2),0, 'normalized')


% apply bandpass
% bandCenters =  10000; % Hz
% band_freq = [-500 500]+bandCenters;
% rir = bandpass(rir,band_freq,fs);
% noise = bandpass(noise,band_freq,fs);

signal = rir + noise;


% estimate noise level
noiseLevel_est = sqrt(mean(noise.^2)) % very similar to ground truth noise level
% noiseLevel = mean(noiseLevel_est);

% compute correlation
r_corr =  movcorr(signal(:,1),signal(:,2),winLen);

e_sig = movsum(signal(:,1).^2,winLen)./winLen;
e_noise = movsum(noise.^2,winLen)./winLen;

% A misestimate of the noise level by a factor can lead to correlation
% between the signal energy and the resulting correlation
% ideal factor is 1
noise_misestimation_factor = 0.9%1.3;
e_rir = e_sig - noiseLevel.^2 * noise_misestimation_factor;
e_rir_true = movsum(rir(:,1).^2,winLen)./winLen;


% for example estimating the noise on a short window can lead to such a
% factor
sqrt(mean(noise(1:1000).^2)) ./ noiseLevel


r_snr = e_rir ./ e_sig;

r_true = (rirLevel/1).^2 ./ ((rirLevel/1).^2 + noiseLevel.^2);

%% robust correlation

med = movmedian(signal,winLen,1);
mad = movmad(signal,winLen,1);

xx = (signal - med) ./ (1 .* mad);
u = xx(:,1) + xx(:,2);
v = xx(:,1) - xx(:,2);

u_med = movmedian(abs(u),winLen,1).^2;
v_med = movmedian(abs(v),winLen,1).^2;
r_med = (u_med - v_med)./(u_med + v_med);
% r_med = r_med / 0.6745;

u_mad = movmad(u,winLen,1).^2;
v_mad = movmad(v,winLen,1).^2;
r_mad = (u_mad - v_mad)./(u_mad + v_mad);
% r_mad = r_mad / 0.6745;

U = sum(signal ./ sqrt(e_sig),2) / sqrt(2);
V = diff(signal ./ sqrt(e_sig),1,2) / sqrt(2);

U_var = movsum(U.^2,winLen)./winLen;
V_var = movsum(V.^2,winLen)./winLen;

r_r = U_var - 1;
r_r = -V_var + 1;

r_r = (U_var - V_var)./(U_var + V_var);



%% plot
figure; hold on; grid on;
plot(r_corr,'LineWidth',1.5)
plot(r_snr * rir_correlation,'LineWidth',1.5)
plot(r_true * rir_correlation,'LineWidth',1.5)
% plot(r_med)
% plot(r_mad)
% plot(r_r)
legend('Signal correlation','SNR-based Correlation','True Correlation')
xlabel('Time (samples)')
ylabel('Normalized Correlation')


figure; hold on; grid on;
plot(r_corr ./ (r_snr * rir_correlation),'LineWidth',1.5)
plot(r_corr ./ (r_true * rir_correlation),'LineWidth',1.5)
legend('Ratio between Signal and SNR-based Correlation','Ratio between Signal and True Correlation')
xlabel('Time (samples)')
ylabel('Normalized Correlation Ratio')

figure; hold on; grid on;
plot(e_rir,'LineWidth',1.5)
plot(rirLevel.^2,'LineWidth',1.5)
legend('Estimated RIR Energy','Specified RIR Energy')
xlabel('Time (samples)')
ylabel('Energy (linear)')

figure; hold on;
bins = linspace(0.8,1.2,1000);
histogram(r_corr ./ r_snr, bins, 'EdgeColor','none')
histogram(r_corr ./ r_true, bins, 'EdgeColor','none')

figure; 
plot([e_rir, e_rir_true, e_sig, e_sig - e_noise(:,1)])
legend('Estimated RIR Energy','True RIR Energy','Signal Energy','Signal-Noise Energy')
xlabel('Time (samples)')
ylabel('Energy (linear)')