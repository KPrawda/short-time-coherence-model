% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-03
clear; clc; close all;
% addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR'
addpath './../Measurements/'

%% Load measurements
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
 
referenceRIR = 5;

switch 1
    case 1
        numRIR = 5;
        filename = 'IR_numClosed_50_numComb_5000_mic_1_sweep_%d.wav';
        direct_delay = 560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay
        minDB = 60;
    case 2
        numRIR = 5;
        filename = 'IR_numClosed_2_numComb_200_mic_1_sweep_%d.wav';
        direct_delay = 560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay
        minDB = 60;
    case 3
        numRIR = 5;
        filename = 'flutter4_setting2_IR_channel_4_sweep_%d.wav';
        direct_delay = 560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        rir = rir(1:3*fs,:);
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay
        minDB = 30;
    case 4
        numRIR = 2;
        filename = 'IR_Synthetic_%d.wav';
        direct_delay = 0;
        fit_onset = 1;
        minDB = 0;
end


rir = rir(fit_onset:end,:); % truncate to time of sound emittance

% rir = rir(:,[2 3]);
% numRIR = 2;

%% additional noise
% rir = rir + randn(size(rir)) * 10^-4 * 3;


%% bandpass filtering and correlation estimation
winLen = 1000;%2^10;
bandCenters =  [8000]; (1:19)*1000; % Hz
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    [rir_band,digitalfilter] = bandpass(rir,band_freq,fs);

    % compute correlation
    rir_band_ref = rir_band(:,referenceRIR);
    rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];

    [cor, energy, snr_temp, r_snr_temp, e_ref] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp;
    snr(:,:,bandIt) = snr_temp;
    meas_energy(:,:,bandIt) = energy;
    meas_cor(:,:,bandIt) = cor;
    meas_energy_ref(:,:,bandIt) = e_ref;

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:,bandIt) = energyDB > min(energyDB) + minDB ;
    
%     mask(round(fs*0.150):end,bandIt)=0;
end
time_cor = (1:size(cor,1)).'/fs; % seconds


%% towards analytical model
[pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs);

%% find divergence
% sinv = @(snr) clip(snr ./ (1 + snr),[eps 1]); % snr to correlation conversion
[divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask, r_snr);

stretched_time = @(div) pred_toastd ./ (sqrt(2)*div);


resampleTime = @(t,x) interp1(t, x, time_cor,'linear','extrap');

for it = 1:numRIR-1
    for bandIt = 1:numel(bandCenters)
        pred_cor_resample(:,it,bandIt) = resampleTime(stretched_time(divergence(it,bandIt)),pred_cor(:,bandIt));
    end
end
%% plot
% plots only first two RIR correlation


figure; hold on; grid on;
offset = 0*(1:size(pred_cor,2));
plot( time_cor, squeeze(meas_cor(:,:,:)),'.')
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample) 
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample .* r_snr,'--','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, r_snr,':','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, mask,'-.')
xlim([0 max(time_cor)])


figure; hold on; grid on;
offset = 0*(1:size(pred_cor,2));
plot( time_cor, squeeze(meas_cor(:,:,:)),'.')
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample .* r_snr,'--','LineWidth',1.5)
xlim([0 max(time_cor)])

figure; hold on; grid on;
plot( time_cor, meas_cor ./ r_snr,'-','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample,'--') 
ylim([0,1])

% figure; hold on; grid on;
% rirInd = (1:size(cor,1))+winLen/2;
% plot( time_cor-winLen/fs/2, db(squeeze(rir_band(rirInd,2:5) - rir_band(rirInd,1))))

% plot( time_cor-winLen/fs/2, db(squeeze(rir_band(rirInd,[1 5]))))
%

jointFactor = meas_cor(:,3) ./ (pred_cor_resample(:,3) .* r_snr(:,3));

figure; hold on; grid on;
offset = 0*(1:size(pred_cor,2));
plot( time_cor, meas_cor ./ sqrt(jointFactor),'.')
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample .* r_snr,'--','LineWidth',1.5)
xlim([0 max(time_cor)])

%% plot
% plot_function_Fig_8_and_9

% figure; hold on; 
% plot(time_cor,rir); 
% set(gca,'ColorOrderIndex',1);
% plot(time_cor,10*(rir(:,[2:5]) - rir(:,1))+0.001);



%% find divergence
function [divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask, snr_cor)
% Fit the divergence 
% Input:
% - pred_toastd = time in seconds of pred_cor
% - pred_cor = predicted correlation of single pulse
% - time_cor = time in seconds of meas_cor
% - meas_cor = measured correlation between measurements
% - mask = high energy region
% - snr_cor = expected correlation based on SNR
[numT, numIR, numBands] = size(meas_cor);
for itIR = 1:numIR
    for itBands = 1:numBands
        cor = meas_cor(:,itIR,itBands) ./ snr_cor(:,itIR,itBands); % SEB: critial change - compensate snr 
        
        m = mask(:,itBands);

        toastd_interp = interp1(pred_cor(:,itBands),pred_toastd,cor,'linear','extrap');

        b = robustfit(time_cor(m),toastd_interp(m),[],[],0);

        divergence(itIR,itBands) = b / sqrt(2); %median(localDivergence(mask(:,itBands)));

    end
end
end


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





%% sliding correlation
function [r_corr, e_sig, snr, r_snr, e_ref] = slidingCorrelation(irRef, ir, winLen)

num_rirs = size(ir,2);

% estimate noise level
noise = ir(round(0.9*end):end,:);
noiseLevel = sqrt(mean(noise.^2));

disp(noiseLevel)

% thresh = mean(noiseLevel)*40;
% irRef(abs(irRef)<thresh) = eps;
% ir(abs(ir)<thresh) = eps;

% compute correlation
% for it = 1:num_rirs
%     r_corr(:,it) =  movcorr(irRef,ir(:,it),winLen);
% 
%     r_cov(:,it) = movsum(irRef.*ir(:,it),winLen)./winLen;
%     e_sig(:,it) = movsum(ir(:,it).^2,winLen)./winLen;
% end

mov = @(x) movsum(x,winLen)./winLen;

win = hann(winLen);
win = rectwin(winLen);
win = blackman(winLen);
win = win / sum(win);
mov = @(x) conv(x,win,'same');

for it = 1:num_rirs
    r_cov(:,it) = mov(irRef.*ir(:,it));
    e_sig(:,it) = mov(ir(:,it).^2);
end
e_ref =  mov(irRef.^2);

% e_sig = movmedian(e_sig,winLen*1,1);
% e_ref = movmedian(e_ref,winLen*1,1);

r_corr = r_cov./sqrt(e_sig.*e_ref);

% r_corr = movmedian(r_corr,winLen*10,1);

% estimate rir energy
e_rir = e_sig - noiseLevel.^2 * 1;

% SNR-based correlation
r_snr = e_rir ./ e_sig;

snr = 1;


med = movmedian(ir,winLen,1);
mad = movmad(ir,winLen,1);

med_ref = movmedian(irRef,winLen,1);
mad_ref = movmad(irRef,winLen,1);


xx = (ir - med) ./ (1 .* mad);
xx_ref = (irRef - med_ref) ./ (1 .* mad_ref);
u = xx_ref + xx;
v = xx_ref - xx;

u_med = movmedian(abs(u),winLen,1).^2;
v_med = movmedian(abs(v),winLen,1).^2;
r_med = (u_med - v_med)./(u_med + v_med);
% r_med = r_med / 0.6745;

u_mad = movmad(u,winLen,1).^2;
v_mad = movmad(v,winLen,1).^2;
r_mad = (u_mad - v_mad)./(u_mad + v_mad);

figure; hold on;
plot(r_corr)
set(gca,'ColorOrderIndex',1);
plot(r_med,'--')
set(gca,'ColorOrderIndex',1);
plot(r_mad,'.-')

ok = 1;

% r_corr = r_med;

% plot correlation ratio
% figure; grid on;
% plot(r_corr ./ r_snr,'LineWidth',1.5)
% ylim([0.99,1])
% xlim([0,20000])
% 
% figure; grid on;
% plot(r_snr,'LineWidth',1.5)
% ylim([0.99,1])
% xlim([0,20000])
% 
% figure; grid on; hold on;
% fs = 44100;
% hp = @(x) highpass(x,150,fs);
% meanremoval = @(x) x - mean(x(6000:8000));
% pp = @(x) meanremoval(hp(x));
% plot(pp(r_snr(:,2)*100),'LineWidth',1.5) 
% plot(pp(r_corr(:,2)),'LineWidth',1.5) 
% % ylim([0.99,1])
% xlim([5000,20000])
% 
% figure; grid on; hold on;
% % plot(r_cov(:,2));
% % plot(e_sig(:,2));
% 
figure; grid on; hold on;
% plot(r_cov(:,:) ./ e_sig(:,:)); 
plot(r_cov(:,:) ./ e_rir(:,:)); 

ylim([0,1.1])

ok = 1;

end