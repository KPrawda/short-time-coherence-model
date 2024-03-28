% Sebastian J. Schlecht, Sunday, 04 December 2022
clear; clc; close all;
% addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR'
addpath './../Measurements/'

%%
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
switch 2
    case 1
        numRIR = 5;
        filename = 'IR_numClosed_50_numComb_5000_mic_1_sweep_%d.wav';
        direct_delay = 0;%560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(rir(:,1));
        fit_onset = 1;%floor(onset) - direct_delay
        minDB = 60;
    case 2
        numRIR = 5;
        filename = 'IR_numClosed_2_numComb_200_mic_1_sweep_%d.wav';
        direct_delay = 560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(rir(:,1));
        fit_onset = floor(onset) - direct_delay
        minDB = 60;
    case 3
        numRIR = 5;
        filename = 'flutter4_setting2_IR_channel_4_sweep_%d.wav';
        direct_delay = 560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(rir(:,1));
        fit_onset = floor(onset) - direct_delay
        minDB = 30;
    case 4
        numRIR = 2;
        filename = 'IR_Synthetic_%d.wav';
        direct_delay = 0;
        fit_onset = 1;
        minDB = 0;
end

for it = 1:numRIR
    [rir(:,it),fs] = audioread(sprintf(filename,it));
end
rir = rir(fit_onset:end,:); % truncate to time of sound emittance

% rir = rir(:,[5 1 2 3 4]);
% rir = rir(1:2*fs,[1 5]); numRIR = 2;

%% bandpass filtering and correlation estimation
winLen = 400;%2^10;
bandCenters =  [10000]; (1:19)*1000; % Hz
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    rir_band = bandpass(rir,band_freq,fs);

    % compute correlation
    [cor, energy, snr_temp, r_snr_temp] = slidingCorrelation2(rir_band, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp(:,2:end);
    snr(:,:,bandIt) = snr_temp(:,2:end);
    meas_energy(:,:,bandIt) = energy(:,2:end);
    meas_cor(:,:,bandIt) = cor(:,2:end);

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:,bandIt) = energyDB > min(energyDB) + minDB ;
end
time_cor = (1:size(cor,1)).'/fs; % seconds




%% towards analytical model
[pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs);

%% find divergence
sinv = @(snr) clip(snr ./ (1 + snr),[eps 1]); % snr to correlation conversion
[divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask, sinv(snr));
% [divergence] = findDivergence2(pred_toastd, pred_cor, time_cor, meas_cor, sinv(snr));

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
plot( time_cor, squeeze(meas_cor(:,:,:)),'LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample) 
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor_resample .* sinv(snr),'--','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, sinv(snr).^100,':','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, mask,'-.')
xlim([0 max(time_cor)])

% figure; hold on; grid on;
% rirInd = (1:size(cor,1))+winLen/2;
% plot( time_cor-winLen/fs/2, db(squeeze(rir_band(rirInd,2:5) - rir_band(rirInd,1))))

% plot( time_cor-winLen/fs/2, db(squeeze(rir_band(rirInd,[1 5]))))
%

%% plot
% plot_function_Fig_8_and_9




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


% %% sliding correlation
% function [xCor, energy, snr] = slidingCorrelation(ir, winLen)
% referenceIR = 1;
% num_rirs = size(ir,2);
% 
% switch 2
%     case 1 % rectangular window
%         win = ones(winLen,1); % / winLen; 
%     case 2 % much smoother than rectangular window
%         win = hann(winLen)*10; 
% end
% % win = win / sum(win);
% 
% for it = 1:num_rirs % for each channel
%     irMean(:,it) = conv(ir(:,it),win/winLen,'valid');
%     xCov(:,it) = conv(ir(:,it) .* ir(:,referenceIR),win,'valid');
%     energy(:,it) = conv(ir(:,it).^2,win,'valid');
%     diffEnergy(:,it) = conv( (ir(:,it) - ir(:,referenceIR)).^2,win,'valid');
% end
% % xCor = xCov ./ sqrt(energy) ./ sqrt(energy(:,referenceIR));
% refMean = irMean(:,referenceIR);
% % SEB: compute actual sample correlation (removing mean)
% xCor = (xCov - winLen .* irMean .* refMean) ./ sqrt(energy - winLen .* irMean.^2) ./ sqrt(energy(:,referenceIR) - winLen .* refMean.^2);
% 
% 
% 
% % SEB: important change -> difference from the background noise was missing
% backgroundNoise = median(energy(round(end*0.9):end,:),[1 2]); % TODO: this background is very hacky and can easily fail
% signalEnergy = mean(energy,2);
% rirEnergy = signalEnergy - backgroundNoise;
% snr = rirEnergy ./ backgroundNoise;
% snr(snr < eps) = eps; % remove negative values
% 
% 
% 
% end



%% sliding correlation
function [xCor, energy, snr, r_snr] = slidingCorrelation2(ir, winLen)
referenceIR = 1;
num_rirs = size(ir,2);

switch 1
    case 1 % rectangular window
        win = ones(winLen,1); % / winLen; 
    case 2 % much smoother than rectangular window
        win = hann(winLen)*10; 
end
% win = win / sum(win);

med = movmedian(ir,winLen,1);
mad = movmad(ir,winLen,1);

xx = ir - med ./ (sqrt(2) .* mad);
u = xx + permute(xx,[1 3 2]);
v = xx - permute(xx,[1 3 2]);

u_med = movmedian(abs(u),winLen,1).^2;
v_med = movmedian(abs(v),winLen,1).^2;
r_med = (u_med - v_med)./(u_med + v_med);

u_mad = movmad(u,winLen,1).^2;
v_mad = movmad(v,winLen,1).^2;
r_mad = (u_mad - v_mad)./(u_mad + v_mad);



for it = 1:num_rirs % for each channel
    irMean(:,it) = conv(ir(:,it),win/winLen,'valid');
    xCov(:,it) = conv(ir(:,it) .* ir(:,referenceIR),win,'valid');
    energy(:,it) = conv(ir(:,it).^2,win,'valid');
    diffEnergy(:,it) = conv( (ir(:,it) - ir(:,referenceIR)).^2,win,'valid');
end
% xCor = xCov ./ sqrt(energy) ./ sqrt(energy(:,referenceIR));
refMean = irMean(:,referenceIR);
% SEB: compute actual sample correlation (removing mean)
xCor = (xCov - winLen .* irMean .* refMean) ./ sqrt(energy - winLen .* irMean.^2) ./ sqrt(energy(:,referenceIR) - winLen .* refMean.^2);


% SEB: important change -> difference from the background noise was missing
backgroundNoise = mean(energy(round(end*0.9):end,:),[1 2]); % TODO: this background is very hacky and can easily fail
signalEnergy = energy; mean(energy,2);
rirEnergy = signalEnergy - backgroundNoise.^2;
rirEnergy(rirEnergy < eps) = eps;
snr = rirEnergy ./ backgroundNoise.^2;
snr(snr < eps) = eps; % remove negative values

r_snr = clip(rirEnergy ./ signalEnergy,[eps,1]);

% SEB: snr removed xcor
% xCor = xCov ./ sqrt(rirEnergy) ./ sqrt(rirEnergy(:,referenceIR));

noiseLevel = mean(ir(round(end*0.7):end,:).^2,[1 2])
e_sig = movsum(ir(:,1).^2,winLen)./winLen;
e_rir = e_sig - noiseLevel;
r_snr = e_rir ./ e_sig;

r_corr = movcorr(ir(:,1),ir(:,2),winLen);
plot(r_corr ./ r_snr,'LineWidth',1.5)
ylim([0,1])
ok = 1;

end