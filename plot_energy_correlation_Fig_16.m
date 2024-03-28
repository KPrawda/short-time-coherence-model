% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-26
clear; clc; close all;
% addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR'
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\Resampled_RIRs\'
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
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
        filename = 'rec1_IR_ch_8_';
        minDB = 30;
        its =  [60, 360, 600];%[60, 360, 600] - 5, 30, 50 min
        numRIR = length(its)+1;
        direct_delay = 560;
        [rir(:, 1), fs] = audioread('rec1_IR_ch_8_1.wav');
        for it = its
            temp = audioread(sprintf('%s%d%s',filename,it, '_resampled.wav'));
            rir(1:length(temp),it) = temp;
        
        end
        rir = rir(:, [1,its]);
        [~,onset] = max(abs(rir(:,1)));
        fit_onset = floor(onset) - direct_delay;
        rir = rir(1:0.9*fs, :);
end


rir = rir(fit_onset:end,:); % truncate to time of sound emittance

%% bandpass filtering and correlation estimation
winLen = 2^10;
bandCenters = (1:19)*1000; % Hz
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    [rir_band,digitalfilter] = bandpass(rir,band_freq,fs);
    
    % compute correlation
    rir_band_ref = rir_band(:,referenceRIR);
    rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];

    [cor, energy, r_snr_temp, e_ref, e_mean] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp;
    meas_energy(:,:,bandIt) = energy;
    mean_energy(:,:,bandIt) = e_mean;
    meas_cor(:,:,bandIt) = cor;
    meas_energy_ref(:,:,bandIt) = e_ref;

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:, bandIt) = energyDB > min(energyDB) + minDB ;
%         mask(round(fs*0.2):end,bandIt)=0;
end
time_cor = (1:size(cor,1)).'/fs; % seconds

%% find divergence
volatility = findVolatility(time_cor, meas_cor.^2, mask, r_snr.^2, bandCenters);
%%

for itIR = 1:numRIR-1
    for bandIt = 1:numel(bandCenters)
        pred_cor_resample(:,itIR,bandIt) = correlationModel(bandCenters(bandIt), time_cor, volatility(itIR,bandIt));       
    end
end
R_2 = pred_cor_resample.*r_snr;

pred_energy= 0.25*(meas_energy_ref + meas_energy + 2*R_2.*sqrt(meas_energy_ref.*meas_energy))./meas_energy_ref;
energy_loss =  mean_energy./meas_energy_ref;
corr_loss = ((1+abs(meas_cor))./2);

energy_lossdB = pow2db(energy_loss);
energy_loss_mod = pow2db(pred_energy);
energy_loss_cor = pow2db(corr_loss);

%% colors
numPlots = numRIR-1;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;

%% plot energy loss (modeled vs measured) for Arni FA over different freq bands
% this one actually does not go into the paper 
f = figure(1); clf; hold on
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band =[10 15 18];
 licz = 1;
for b = band 
     for i =2        
        plot( 1000* time_cor, (squeeze(energy_lossdB(:,i,b))),'-', 'color', cMap2(:,licz), 'LineWidth',1.2)
        plot( 1000* time_cor,(squeeze(energy_loss_mod(:,i,b))) ,'--', 'color', cMap2(:,licz), 'LineWidth',3.5, 'HandleVisibility','off')
     end    
     licz = licz+1 ;
end
ylim([-4.2 0.3])
xlim([0 320])
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Energy loss (dB)', 'Interpreter','latex', 'FontSize',12)
set(gca, 'XTick', 0:50:450, 'fontsize', 12)
lgd  = legend('10 kHz', '15 kHz', '18 kHz',  'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 1);
% lgd.Title.String = [ {'Center'}, { 'frequency (kHz)'}];%$j-i$';
box on
grid on
f.Position(end) = 320;
%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_FA_energy_loss_freq','-dpdf','-r0')

%% plot energy loss (modeled vs measured) for Arni FA over different time separations
% this one goes into the paper
f = figure(2); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band =18;
for b = band 
     for i = 1: numRIR-1
        plot( 1000* time_cor, (squeeze(energy_lossdB(:,i,b))),'-', 'color', cMap2(:,i), 'LineWidth',1.2)
        plot( 1000* time_cor,(squeeze(energy_loss_mod(:,i,b))) ,'--', 'color', cMap2(:,i), 'LineWidth',3.5, 'HandleVisibility','off')       
     end    
end
ylim([-4.2 0.3])
xlim([0 270])
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Energy loss (dB)', 'Interpreter','latex', 'FontSize',12)
set(gca, 'XTick', 0:50:300, 'fontsize', 12)
lgd  = legend('5 min', '30 min', '50 min',   'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 3);
lgd.Title.String = [ {'Time between measurements'}];%$j-i$';
box on
grid on
f.Position(end) = 320;

%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_FA_energy_loss_time','-dpdf','-r0')
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
function [r_corr, e_sig, r_snr, e_ref, e_mean] = slidingCorrelation(irRef, ir, winLen)

num_rirs = size(ir,2);

% estimate noise level
noise = irRef(round(0.9*end):end,:);
noiseLevel = sqrt(mean(noise.^2,1));
size(noiseLevel)
% disp(noiseLevel)

% compute correlation
mov = @(x) movsum(x,winLen)./winLen;

% win = hann(winLen);
win = rectwin(winLen);
% win = blackman(winLen);
win = win / sum(win);
mov = @(x) conv(x,win,'same');

% avgEnergy = abs(avg_s).^2;
% avgEnergy = movsum(avgEnergy,freqWin,1);
% avgEnergy = movsum(avgEnergy,timeWin,2);

for it = 1:num_rirs

    r_cov(:,it) = mov(irRef.*ir(:,it));
    e_sig(:,it) = mov(ir(:,it).^2);

    sig = [irRef, ir(:,it)];
    avg_s = mean(sig,2);
    avgEnergy = abs(avg_s).^2;
    avgEnergy = mov(avgEnergy);
    e_mean(:, it) =  avgEnergy;
end
e_ref =  mov(irRef.^2);
r_corr = r_cov./sqrt(e_sig.*e_ref);

% estimate rir energy
e_rir = max(e_sig - noiseLevel.^2 * 1,0);

% SNR-based correlation
r_snr = e_rir ./ e_sig;


end