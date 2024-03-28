%% plot figures 8, 9 and 11 from the manuscript
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024

% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-26
% updated by K. Prawda 12.02.2024
% PLOT ONLY THE CORRELATION FOR ARNI 1 DATASET

%% housekeeping
clear; clc; close all;
addpath './RIRs/'
%% Load measurements
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
 
referenceRIR = 1;

switch 1
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

    [cor, energy, r_snr_temp, e_ref] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp;
    meas_energy(:,:,bandIt) = energy;
    meas_cor(:,:,bandIt) = cor;
    meas_energy_ref(:,:,bandIt) = e_ref;

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:, bandIt) = energyDB > min(energyDB) + minDB ;
    %     mask(round(fs*0.150):end,bandIt)=0;
end
time_cor = (1:size(cor,1)).'/fs; % seconds

%% find divergence
volatility = findVolatility(time_cor, meas_cor, mask, r_snr, bandCenters);
%% coherence model
for itIR = 1:numRIR-1
    for bandIt = 1:numel(bandCenters)
        pred_cor_resample(:,itIR,bandIt) = correlationModel(bandCenters(bandIt), time_cor, volatility(itIR,bandIt));
    end
end
%% colors
numPlots = numRIR-1;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;
%% plot Arni 1 - measured, modeled, expected and so on, figure 8

f=figure(1); clf; hold on
subplot(2,1,1); hold on
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band = 10; % which band to plot
nIR = 4; % which RIR to plot
%measured correlation
plot(1000*time_cor,meas_cor(:,nIR,band).^2 ,'-', 'color', cMap2(:, 1), 'LineWidth',3)
% expected correlation from the noise
plot(1000*time_cor, squeeze(r_snr(:, nIR, band)).^2, 'color', cMap2(:, 2), 'LineWidth',2)
%model prediction
plot( 1000*time_cor, pred_cor_resample(:,nIR,  band).^2 ,'-.',  'color', cMap2(:, 3), 'LineWidth',2)
% plot exptected + prediction
plot(  1000*time_cor,(r_snr(:, nIR, band).*pred_cor_resample(:,nIR, band)).^2,'--' ,'color', cMap2(:, 4), 'LineWidth',2)
x = correlationModel(10000,time_cor,volatility(4,1)).*r_snr(:, 4,1);
% plot(1000*time_cor,x, 'k');
plot(time_cor*1000, 2*mask(:, band),'--','color', 'k' ,'LineWidth',0.5, 'HandleVisibility','off'); hold on

ylim([0.9 1.01])
xlim(1000*[0 0.5]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Coherence', 'Interpreter','latex', 'FontSize',12)
set(gca, 'ytick', .9:.02:1, 'xtick', 0:100:600, 'FontSize',12)
lgd  = legend('$\widehat{\rho}^{~2}_{x, x^\prime}$', '$\widehat{\rho}^{~2}_{\textup{SNR}}$', '$R^2_{\vartheta}(t, f)$', '$\widehat{\rho}^{~2}_{\textup{SNR}} \times R^2_{\vartheta}(t, f)$', 'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
box on
grid on
%% plot SNR
subplot(2,1,2);  hold on
b = 10;
nIR = 1;

band_freq = [-500 500]+bandCenters(b);
[rir_band,digitalfilter] = bandpass(rir,band_freq,fs);

% compute correlation
rir_band_ref = rir_band(:,referenceRIR);
rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];

[~, energyb, ~, ~] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);

enoise = median(energyb(0.9*length(energyb):end, nIR));
snr = (energyb(:, nIR)-enoise)./enoise;

plot(time_cor*1000, 10*log10(abs(snr)),'-','color', 'k' ,'LineWidth',1.5); hold on
plot(time_cor*1000, 100*mask(:, b)-20,'--','color', 'k' ,'LineWidth',0.5, 'HandleVisibility','off'); hold on

xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('SNR (dB)', 'Interpreter','latex', 'FontSize',12)

% xline(27)
xlim([0 500])
ylim([-5  75])
grid on
box on
set(gca, 'FontSize', 12)
% legend('$\widehat{\rho}_{x, x^prime}$', 'SNR', 'Interpreter','latex', 'fontsize', 12, 'location', 'southwest', 'numcolumns',3)
%% print figure 1
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_different_correlations','-dpdf','-r0')
%% plot Arni - different bands, modeled and measured coherence, figure 9
f = figure(2); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band =[10 19];
 licz = 1;
for b = band   
    plot( 1000*time_cor, squeeze(meas_cor(:,nIR,b)).^2 ,'-', 'color', cMap2(:,licz), 'LineWidth',1.2)
    plot(  1000*time_cor, (squeeze(r_snr(:, nIR, b)).*pred_cor_resample(:,nIR, b)).^2 ,'--' ,'color', cMap2(:, licz), 'LineWidth',3)
    licz = licz+3;
end

ylim([-0.3 1.03])
xlim(1000*[0 0.65]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Coherence', 'Interpreter','latex', 'FontSize',12)
lgd  = legend('$\widehat{\rho}^{~2}_{x, x^\prime}$, 10 kHz', '$\widehat{\rho}^{~2}_{\textup{SNR}} \times R^2_{\vartheta}$, 10 kHz', '$\widehat{\rho}^{~2}_{x, x^\prime}$, 19 kHz', '$\widehat{\rho}^{~2}_{\textup{SNR}} \times R^2_{\vartheta}$, 19 kHz','location', 'south', 'interpreter', 'latex', 'fontsize', 11, 'numcolumns', 2);
box on
grid on
set(gca, 'FontSize', 12)
f.Position(end) = 320;
%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_correlation_Jan2024','-dpdf','-r0')

%% get the correlation coefficient between the calculated correlation curves and the fitted curves
cor_PCC = zeros(2,2,numRIR-1,numel(bandCenters));
for i = 1:numRIR-1
    for it = 1: numel(bandCenters)
        firstI = find(mask(:, it) == 1, 1, 'first');
        corL = find(mask(firstI+1:end, it) == 0, 1, 'first');
        nc_len = corL;
        corr_temp =  squeeze(r_snr(:, i, it)).*pred_cor_resample(:,i, it);
       
        cor_PCC (:, :, i, it) = corrcoef([meas_cor(firstI:firstI+nc_len-1,i, it), corr_temp(firstI:firstI+nc_len-1)]);
    end
end
%% plot the model-signal correlation coefficient for the fits, figure 11
f = figure(3); clf; hold on

for i =1:numRIR-1 
    plot(bandCenters(1:19), squeeze(cor_PCC(1,2,i, :)), '-.o', 'color', cMap2(:, i), 'LineWidth',1, 'MarkerFaceColor',cMap2(:, i))
end 

ylim([0.82 1])
xlim(1000*[0.5 19.5]);
xlabel('Center frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('Model-signal correlation', 'Interpreter','latex', 'FontSize',12)
set(gca,  'xtick', 1000*[1:2:20], 'xticklabel', 1:2:20,  'FontSize',12)
box on
f.Position(end) = 250;
lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 4);
lgd.Title.String = [ {'Time between measurements'}];
%% print the correlation coefficient figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'corr_coef_Arni','-dpdf','-r0')
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