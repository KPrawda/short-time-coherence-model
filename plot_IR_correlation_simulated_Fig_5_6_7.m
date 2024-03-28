%% plot figures 5, 6, and 7 from the manuscript
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024

% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-26
% updated by K. Prawda 12.02.2024
% PLOT ONLY THE SIMULATION RESULTS, FIGS 5 AND 6
%% houskeeping
clear; clc; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%% Load measurements
referenceRIR = 1;

numRIR = 2;
filename = 'IR_Synthetic_%d.wav';
for itIR = 1:numRIR
    [rir(:,itIR),fs] = audioread(sprintf(filename,itIR));
end
direct_delay = 0;
fit_onset = 1;
minDB = 0;

rir = rir(fit_onset:end,:); % truncate to time of sound emittance
%% bandpass filtering and correlation estimation
winLen = 2^10;
bandCenters = (1:20)*1000; % Hz
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    [rir_band,digitalfilter] = bandpass(rir,band_freq,fs);

    % compute correlation
    rir_band_ref = rir_band(:,referenceRIR);
    rir_band_other = rir_band; rir_band_other(:,referenceRIR) = [];

    [cor, energy, r_snr_temp, e_ref] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp.^2;
    meas_energy(:,:,bandIt) = energy;
    meas_cor(:,:,bandIt) = cor.^2;
    meas_energy_ref(:,:,bandIt) = e_ref;

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:, bandIt) = energyDB > min(energyDB) + minDB ;
    %     mask(round(fs*0.150):end,bandIt)=0;
end
time_cor = (1:size(cor,1)).'/fs; % seconds

%% find divergence
volatility = findVolatility(time_cor, sqrt(meas_cor), mask, r_snr, bandCenters);
%% apply coherence model
for itIR = 1:numRIR-1
    for bandIt = 1:numel(bandCenters)
        pred_cor_resample(:,itIR,bandIt) = correlationModel(bandCenters(bandIt), time_cor, volatility(itIR,bandIt));
    end
end
%% colors
numPlots =3;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;
%% plot simulated over time, figure 5
nIR = 1;
licz = 1;
f = figure(2); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band =[ 5 10 20];
for b = band 
    plot( 1000*time_cor, squeeze(meas_cor(:,nIR,b)) ,'-', 'color', cMap2(:,licz), 'LineWidth',1.5)    
    plot(  1000*time_cor, pred_cor_resample(:,nIR, b).^2,'--' ,'color', cMap2(:, licz), 'LineWidth',2, 'HandleVisibility','off')
    licz = licz+1;
end

set(gca, 'fontsize', 12)
% ylim([-0.05 1.03])
ylim([0.55 1.01])
xlim(1000*[0 1]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('$R^2_\vartheta (\tau, f)$', 'Interpreter','latex')
lgd  = legend('5 kHz', '10 kHz','20 kHz','location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 1);
box on
grid on
%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Simulated_corr_over_time','-dpdf','-r0')

%% plot simulated over frequencies, figure 6
nIR = 1;
licz = 1;
f = figure(1); clf; hold on
t = [0.25 0.5 1]*fs;

for t_ = t
    plot( bandCenters./1000, squeeze(meas_cor(t_,nIR,:)) ,'-', 'color', cMap2(:,licz),  'LineWidth',1.5)
    plot(  bandCenters./1000, squeeze(pred_cor_resample(t_,nIR, :)).^2,'--', 'color', cMap2(:,licz), 'LineWidth',2, 'HandleVisibility','off')
    licz = licz+1;
end


set(gca, 'fontsize', 12)
% ylim([-0.05 1.03])
xlim([1 20]);
xlabel('Frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('$R^2_\vartheta (\tau, f)$', 'Interpreter','latex')
lgd  = legend('250 ms', '500 ms','1000 ms','location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 1);
box on
grid on
% f.Position(end) = 320;

%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Simulated_corr_over_frequency','-dpdf','-r0')

%% get the correlation coefficient between the calculated correlation curves and the fitted curves
cor_PCC = zeros(2,2,numRIR-1,numel(bandCenters));
for i = 1:numRIR-1
    for it = 1: numel(bandCenters)
               
        cor_PCC (:, :, i, it) = corrcoef([meas_cor(:,i, it), pred_cor_resample(:,i, it).^2]);
    end
end
%% plot the model-signal correlation coefficient for the fits
f = figure(3); clf; hold on

for i =1:numRIR-1 
    plot(bandCenters, squeeze(cor_PCC(1,2,i, :)), '-.o', 'color', cMap2(:, i), 'LineWidth',1, 'MarkerFaceColor',cMap2(:, i))
end 

ylim([0.75 1])
xlim(1000*[0.5 19.5]);
xlabel('Center frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('Model-signal correlation', 'Interpreter','latex', 'FontSize',12)
set(gca,  'xtick', 1000*[1:2:20], 'xticklabel', 1:2:20,  'FontSize',12)
box on
f.Position(end) = 250;
lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'south', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 4);
lgd.Title.String = [ {'Time between measurements'}];

%% plot figure 7 - measured volatility vs target
f = figure(4); clf; hold on

plot(bandCenters./1000, volatility, 'LineWidth',2)

xlabel('Frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('$\vartheta (s/\sqrt{s})$', 'Interpreter','latex', 'FontSize',12)

yline(5*10^-6, '--','LineWidth',2)
xlim([1 20])
% ylim([2.75 3.08]*10^-5)
grid on
box on
set(gca, 'FontSize', 12, 'xtick', [1 2:2:20])
legend('Estimated', 'Ground truth', 'Interpreter','latex', 'fontsize', 12, 'location', 'northeast', 'numcolumns',3)

f.Position(end) = 270;
%% print the correlation vs SNR figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'volatility_simulated','-dpdf','-r0')
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
        snr = ones(size(snr_cor(m,itIR,itBands)));
        F = fb(itBands);

        % l1 loss
        loss_fun = @(volatility) sum(abs(correlationModel(F,T,exp(volatility)).*snr - cor));

        % do the fitting in log(volatility) for better scaling
        volatilityMin = -50;
        volatilityMax = -3;
        options = optimset('TolX',0.01);
        %%
        vol = exp(fminbnd(loss_fun,volatilityMin,volatilityMax,options));
        %%
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