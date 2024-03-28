% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-26
clear; clc; close all;
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR'
% addpath './../Measurements/'

%% Load measurements
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
 
referenceRIR = 1;

switch 2
    case 1
        numRIR = 5;
        filename = 'IR_numClosed_50_numComb_5000_mic_1_sweep_%d.wav';
        direct_delay = 560;
        for itIR = 1:numRIR
            [rir(:,itIR),fs] = audioread(sprintf(filename,itIR));
        end
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay;
        minDB = 60;
    case 2
        numRIR = 5;
        filename = 'IR_numClosed_2_numComb_200_mic_1_sweep_%d.wav';
        direct_delay = 560;
        for itIR = 1:numRIR
            [rir(:,itIR),fs] = audioread(sprintf(filename,itIR));
        end
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay;
        minDB = 60;
    case 3
        numRIR = 5;
        filename = 'flutter4_setting2_IR_channel_4_sweep_%d.wav';
        direct_delay = 560;
        for itIR = 1:numRIR
            [rir(:,itIR),fs] = audioread(sprintf(filename,itIR));
        end
        rir = rir(1:3*fs,:);
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay;
        minDB = 30;
    case 4
        numRIR = 2;
        filename = 'IR_Synthetic_%d.wav';
        direct_delay = 0;
        fit_onset = 1;
        minDB = 0;
end


rir = rir(fit_onset:end,:); % truncate to time of sound emittance

%% bandpass filtering and correlation estimation
winLen = 2^10;
bandCenters = 5000; (2:5:19)*1000; % Hz
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    [rir_band,digitalfilter] = bandpass(rir,band_freq,fs);

    % compute correlation
    rir_band_ref = rir_band(:,referenceRIR);
    rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];

    [cor, energy, r_snr_temp, e_ref] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp;
    meas_energy(:,:,bandIt) = energy;
    meas_cor(:,:,bandIt) = cor;
    meas_energy_ref(:,:,bandIt) = e_ref;

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:,bandIt) = energyDB > min(energyDB) + minDB ;
    %     mask(round(fs*0.150):end,bandIt)=0;
end
time_cor = (1:size(cor,1)).'/fs; % seconds

%% find divergence
divergence = findDivergence(time_cor, meas_cor, mask, r_snr, bandCenters);

divergence * 1000000 % unit = microseconds per second

for itIR = 1:numRIR-1
    for bandIt = 1:numel(bandCenters)
        pred_cor(:,itIR,bandIt) = correlationModel(bandCenters(bandIt), time_cor, divergence(itIR,bandIt));
    end
end




%% plot
figure; hold on; grid on;
plot( time_cor, squeeze(meas_cor(:,:,:)),'.')
set(gca,'ColorOrderIndex',1);
% plot( time_cor, pred_cor_resample) 
set(gca,'ColorOrderIndex',1);
plot( time_cor, pred_cor .* r_snr,'--','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, r_snr,':','LineWidth',1.5)
set(gca,'ColorOrderIndex',1);
plot( time_cor, mask,'-.')
xlim([0 max(time_cor)])

%% plot
% plot_function_Fig_8_and_9


%% find divergence
function [divergence] = findDivergence(time_cor, meas_cor, mask, snr_cor, fb)
% Fit the divergence 
% Input:
% - time_cor = time in seconds of meas_cor
% - meas_cor = measured correlation between measurements
% - mask = high energy region
% - snr_cor = expected correlation based on SNR
% - fb = center frequencies
% - fs = sampling frequencies

[numT, numIR, numBands] = size(meas_cor);
for itIR = 1:numIR
    for itBands = 1:numBands
        m = mask(:,itBands);

        T = time_cor(m);
        cor = meas_cor(m,itIR,itBands);
        snr = snr_cor(m,itIR,itBands);
        F = fb(numBands);

        % l1 loss
        loss_fun = @(divergence) sum(abs(correlationModel(F,T,exp(divergence)).*snr - cor));

        % do the fitting in log(divergence) for better scaling
        divergenceMin = -50;
        divergenceMax = -3;
        options = optimset('TolX',0.01);
        div = exp(fminbnd(loss_fun,divergenceMin,divergenceMax,options));

        divergence(itIR,itBands) = div;

        figure; hold on;
        plot(T,cor)
        plot(T,correlationModel(F,T,div).*snr);
        
        ok = 1;
    end
end
end

%% correlation model
function [pred_cor,pred_toastd] = correlationModel(F,T,divergence)
    magicFactor = 20;
    pred_cor = exp( - magicFactor * (F .* sqrt(T)*divergence).^2 );
    pred_toastd = T*divergence;
end


%% sliding correlation
function [r_corr, e_sig, r_snr, e_ref] = slidingCorrelation(irRef, ir, winLen)

num_rirs = size(ir,2);

% estimate noise level
noise = ir(round(0.9*end):end,:);
noiseLevel = sqrt(mean(noise.^2));

disp(noiseLevel)

% compute correlation
mov = @(x) movsum(x,winLen)./winLen;

% win = hann(winLen);
win = rectwin(winLen);
% win = blackman(winLen);
win = win / sum(win);
mov = @(x) conv(x,win,'same');

for it = 1:num_rirs
    r_cov(:,it) = mov(irRef.*ir(:,it));
    e_sig(:,it) = mov(ir(:,it).^2);
end
e_ref =  mov(irRef.^2);
r_corr = r_cov./sqrt(e_sig.*e_ref);

% estimate rir energy
e_rir = e_sig - noiseLevel.^2 * 1;

% SNR-based correlation
r_snr = e_rir ./ e_sig;


end