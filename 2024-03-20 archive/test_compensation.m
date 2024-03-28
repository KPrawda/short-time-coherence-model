%% compensate for the time-variance-induced loss of energy

clear; clc; close all;
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR'
%%
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
switch 3
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
        fit_onset = floor(onset) - direct_delay;
        minDB = 60;
    case 3
        numRIR = 5;
        filename = 'flutter4_setting2_IR_channel_4_sweep_%d.wav';
        direct_delay = 560;
        for it = 1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(rir(:,1));
        fit_onset = floor(onset) - direct_delay;
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
% rir = circshift(rir,-1,2); % reorder the sweeps​
% ​
%% bandpass filtering and correlation estimation
winLen = 1000;%2^10;
bandCenters =  (1:19)*1000; % Hz
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    rir_band = bandpass(rir,band_freq,fs);

    % compute correlation
    [cor, energy] = slidingCorrelation(rir_band, winLen);

    meas_cor(:,:,bandIt) = cor(:,2:end);

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:,bandIt) = energyDB > min(energyDB) + minDB ;
    %     mask(1:(winLen/2 + direct_delay),:) = 0;
end
time_cor = (1:size(cor,1)).'/fs; % seconds

%% towards analytical model
[pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs);

%% find divergence
[divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask);

%% get the correlation coefficient between the calculated correlation curves and the fitted curves
clear new_cor

for i = 1:4
    for it = 1: bandIt
        firstI = find(mask(:, it) == 1, 1, 'first');
        firstI = firstI+500;
        corL = find(mask(firstI+1:end, it) == 0, 1, 'first');

        pred_time(:, i, it) = pred_toastd ./ (sqrt(2)*divergence(i,it));
        [~, loc] = min(abs(pred_time(1:end, i, it) - corL/fs));
        new_cor = interp1(pred_time(1:loc, i, it), pred_cor(1:loc,  it), time_cor(firstI:corL));
%         new_cor = new_cor(~isnan(new_cor));
%         nc_len = length(new_cor);
       
    end
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
function [xCor, energy] = slidingCorrelation(ir, winLen)
referenceIR = 1;
num_rirs = size(ir,2);
% win = ones(frameLen,1) / (frameLen);
win = hann(winLen); % much smoother than rectangular window
win = win / sum(win);

for it = 1:num_rirs % for each channel
    xCov(:,it) = conv(ir(:,it) .* ir(:,referenceIR),win,'valid');
    energy(:,it) = conv(ir(:,it).^2,win,'valid');
    diffEnergy(:,it) = conv( (ir(:,it) - ir(:,referenceIR)).^2,win,'valid');
end
xCor = xCov ./ sqrt(energy) ./ sqrt(energy(:,referenceIR));

end