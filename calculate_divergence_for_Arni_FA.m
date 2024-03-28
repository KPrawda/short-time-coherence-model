% estimate the divergence for Arni RIRs
clear; clc; close all;

% %%
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\Resampled_RIRs\';
% addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Forum Acusticum 2023\RIRs\';
%% read the corresponding IR from the Arni database

folder='C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\Resampled_RIRs\';
% folder='C:\Users\prawdak1\Dropbox (Aalto)\Projects\Forum Acusticum 2023\RIRs\';
audio_files=dir(fullfile(folder,'*.wav'));

file = cell(numel(audio_files),1);
for k=1:numel(audio_files)
  filename=audio_files(k).name;
  if str2double(filename(4)) == 1 || str2double(filename(4)) == 2
      file{k} = filename;    
  end

end
% file = file';
%%
winLen = 2^10;
bandCenters =  (1:19)*1000; % Hz
minDB = 60;
%%
referenceRIR = 1;
[ref_rir, fs]  = audioread(file{1}); 
ref_rir(1.5*fs+1:end, :) = [];

%
for n = 1001:1293%10%904:1291
    idx = [];
    rir(:, 1) = ref_rir;
    clear meas_cor mask
    if n < 632
        filN = sprintf('%s%d%s', 'rec1_IR_ch_8_', n, '_resampled.wav');
    elseif n > 631
        filN = sprintf('%s%d%s', 'rec2_IR_ch_8_', n-631, '_resampled.wav');
    end
    for k = 1:numel(file)
        if contains(file{k}, filN)
            idx = k;
            break;
        end
    end
    if isempty(idx)
        continue
    end
    direct_delay = 0;
    rir_temp  = audioread(file{idx}); 
    len = min(length(rir_temp), 1.5*fs);
    rir(1:len, 2)=rir_temp(1:len);
    %%
    [value, location] = max((rir));
    if diff(location) > 0
        rir(:, 2) = [rir(diff(location)+1:end, 2); rir(end -diff(location)+1:end, 2)];
    elseif diff(location) < 0
        rir(:, 1) = [rir(-diff(location)-2:end, 1); rir(end +diff(location)+4:end, 1)];
    end


    for bandIt = 1:numel(bandCenters)
        band_freq = [-500 500]+bandCenters(bandIt);
        [rir_band,digitalfilter] = bandpass(rir,band_freq,fs);
    
        % compute correlation
        rir_band_ref = rir_band(:,referenceRIR);
        rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];
    
        [cor, energy, r_snr_temp, e_ref] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
        
        r_snr(:, bandIt) = r_snr_temp;
        meas_energy(:,bandIt) = energy;
        meas_cor(:,bandIt) = cor;
        meas_energy_ref(:,bandIt) = e_ref;
    
        % find valid part
%         energyDB = db(mean(energy,2));
        energyDB = db(energy);
        mask(:, bandIt) = energyDB > min(energyDB) + minDB ;
        mask(round(fs*0.2):end,bandIt)=0;
    end
    time_cor = (1:size(cor,1)).'/fs; % seconds
    %% find divergence
    volatility(n,:) = findVolatility(time_cor, meas_cor, mask, r_snr, bandCenters);
    %% model coherence
    for bandIt = 1:numel(bandCenters)
       pred_cor_resample(:,bandIt) = correlationModel(bandCenters(bandIt), time_cor, volatility(n,bandIt));
    end   
end
%% save the results to a new variable

% save('Arni_div_new.mat', 'v_s', 'numComb');
% load('Arni_div_new.mat');
save('Arni_volatility_FA.mat', 'volatility');
% load('Arni_div_new.mat');

%% plot test
band = 10
figure(1); clf; hold on
plot(time_cor, meas_cor(:, band))
plot(time_cor, mask(:, band))
corr_temp = interp1( pred_toastd ./ (sqrt(2)*v_s(band, n)), pred_cor(:, band) , time_cor);
plot(time_cor, corr_temp)
plot(time_cor, exp_corr(:, band,1).*corr_temp)
%% find divergence
function [volatility] = findVolatility(time_cor, meas_cor, mask, snr_cor, fb)
% Fit the volatility 
% Input:
% - time_cor = time in seconds of meas_cor
% - meas_cor = measured correlation between measurements
% - mask = high energy region
% - snr_cor = expected correlation based on SNR
% - fb = center frequencies
% - fs = sampling frequencies

[numT, numBands] = size(meas_cor);

    for itBands = 1:numBands
        m = mask(:,itBands);

        T = time_cor(m);
        cor = meas_cor(m,itBands);
        snr = snr_cor(m,itBands);
        F = fb(itBands);

        % l1 loss
        loss_fun = @(volatility) sum(abs(correlationModel(F,T,exp(volatility)).*snr - cor));

        % do the fitting in log(volatility) for better scaling
        volatilityMin = -50;
        volatilityMax = -3;
        options = optimset('TolX',0.01);
        vol = exp(fminbnd(loss_fun,volatilityMin,volatilityMax,options));

        volatility(itBands) = vol;

%         figure; hold on;
%         plot(T,cor)
%         plot(T,correlationModel(F,T,vol).*snr);
%         plot(T, m(m))
        
        ok = 1;
    end

end

%% correlation model
function [pred_cor,pred_toastd] = correlationModel(F,T,volatility)
    magicFactor = 20;
    pred_cor = exp( - magicFactor * (F .* sqrt(T)*volatility).^2 );
    pred_toastd = T*volatility;
end


%% sliding correlation
function [r_corr, e_sig, r_snr, e_ref] = slidingCorrelation(irRef, ir, winLen)

num_rirs = size(ir,2);

% estimate noise level
noise = irRef(round(0.9*end):end,:);
noiseLevel = sqrt(mean(noise.^2));

% disp(noiseLevel)

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