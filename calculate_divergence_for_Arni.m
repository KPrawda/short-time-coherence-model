% estimate the divergence for Arni RIRs
clear; clc; close all;

%%
nTrials = 150;
rng(1)
numComb = floor(rand(1, nTrials)*5341)+1;
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Mesurements Arni\IR_Arni_upload\IR_Arni_upload';
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Mesurements Arni';

%% get the total number of files
dat = readmatrix('combinations4.csv');
combinations = dat(2:end, :); %first row is the number of panel

[~, numPanels]  = size(combinations);
numClosed = 55-sum(combinations, 2);

%% read the corresponding IR from the Arni database
folder='C:\Users\prawdak1\Dropbox (Aalto)\Mesurements Arni\IR_Arni_upload\IR_Arni_upload';
audio_files=dir(fullfile(folder,'*.wav'));

file = cell(numel(audio_files),1);
for k=1:numel(audio_files)
  filename=audio_files(k).name;

      file{k} = filename;     

end
file = file';
fs = 44100;
%%
winLen = 2^10;
bandCenters =  (1:19)*1000; % Hz
minDB = 60;
volatility = zeros(nTrials, 5,4,numel(bandCenters));
%%
for n =11:150%:nTrials
    rirInd = 1:5;
    clear rir_band meas_cor mask
    filN = sprintf('%s%d%s%d%s%d%s%d%s',"IR_numClosed_",numClosed(numComb(n)+1),"_numComb_", numComb(n),"_mic_1");
    direct_delay = 560;
    rir = zeros(2.4*fs, 5);
    idx = find(contains(file,filN));
        for i = 1:length(idx)
            rir(:, i) = audioread(file{idx(i)}); 
        end
%% Get starting point of curve fit
        [~,onset] = max(rir(:,1));
        fit_onset = floor(onset) - direct_delay;
        rir = rir(fit_onset:end,:);
    % divergence_start = 0;
    isValidRIR = sum(rir,1);
    rirInd(isValidRIR == 0) = [];
    %% per band fit
    for referenceRIR = rirInd
        for bandIt = 1:numel(bandCenters)
            band_freq = [-500 500]+bandCenters(bandIt);
            [rir_band,digitalfilter] = bandpass(rir,band_freq,fs);
         
    
    %%
    
        % compute correlation
        rir_band_ref = rir_band(:,referenceRIR);
        rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];
        %%
            [cor, energy, r_snr_temp, e_ref] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
        %%
            r_snr(1:length(r_snr_temp),:,bandIt) = r_snr_temp;
            meas_energy(1:length(r_snr_temp),:,bandIt) = energy;
            meas_cor(1:length(r_snr_temp),:,bandIt) = cor;
            meas_energy_ref(1:length(r_snr_temp),:,bandIt) = e_ref;
                    % find valid part
            energyDB = db(mean(energy,2));
            % energyDB = db(energy);
            mask(:, bandIt) = energyDB > min(energyDB) + minDB ;
            % mask(round(fs*0.2):end,bandIt)=0;
            time_cor = (1:size(cor,1)).'/fs; % seconds
        end
        %% find divergence
            volatility(n,referenceRIR, :, :)  = findVolatility(time_cor, meas_cor, mask, r_snr, bandCenters);
            
    end 
    
        
end
%%
save('Arni_volatility.mat', 'volatility', 'numComb');
% load('Arni_div_new.mat');


%% plot
figure; hold on; grid on;
plot(bandCenters,divergence')
xlabel('Frequency Band (Hz)');
ylabel('Divergence (seconds)');

% plots only first two RIR correlation
figure; hold on; grid on;
offset = 0*(1:size(pred_cor,2));
plot( time_cor, squeeze(meas_cor(:,1,:)) + offset)
set(gca,'ColorOrderIndex',1);
plot( pred_toastd ./ (sqrt(2)*divergence(1,:)), pred_cor + offset)
xlim([0 max(time_cor)])


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

[numT, numIR, numBands] = size(meas_cor);
for itIR = 1:numIR
    for itBands = 1:numBands
        m = mask(:,itBands);

        T = time_cor(m);
        cor = meas_cor(m,itIR,itBands);
        snr = snr_cor(m,itIR,itBands);
        F = fb(itBands);

        % l1 loss
        loss_fun = @(volatility) sum(abs(correlationModel(F,T,exp(volatility)).*snr - cor));

        % do the fitting in log(volatility) for better scaling
        volatilityMin = -50;
        volatilityMax = -3;
        options = optimset('TolX',0.01);
        vol = exp(fminbnd(loss_fun,volatilityMin,volatilityMax,options));

        volatility(itIR,itBands) = vol;

%         figure; hold on;
%         plot(T,cor)
%         plot(T,correlationModel(F,T,vol).*snr);
%         plot(T, m(m))
%         
        ok = 1;
    end
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