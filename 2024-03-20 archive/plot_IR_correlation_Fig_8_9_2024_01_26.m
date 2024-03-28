% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-26
clear; clc; %close all;
% addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR'
addpath './RIRs/'

%% Load measurements
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
 
referenceRIR = 1;

switch 4
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
        for itIR = 1:numRIR
            [rir(:,itIR),fs] = audioread(sprintf(filename,itIR));
        end
        direct_delay = 0;
        fit_onset = 1;
        minDB = 0;
end


rir = rir(fit_onset:end,:); % truncate to time of sound emittance

%% bandpass filtering and correlation estimation
winLen = 2^10;
bandCenters = 10000; (2:5:19)*1000; % Hz
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
volatility = findVolatility(time_cor, meas_cor, mask, r_snr, bandCenters);

for itIR = 1:numRIR-1
    for bandIt = 1:numel(bandCenters)
        pred_cor_resample(:,itIR,bandIt) = correlationModel(bandCenters(bandIt), time_cor, volatility(itIR,bandIt));
    end
end




%% plot
figure; hold on; grid on;
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

%% colors
numPlots = 4%numRIR-1;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;
%% plot Arni - measured, modeled, expected and so on 

f=figure(3); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band = 1; % which band to plot
nIR = 1; % which RIR to plot
%measured correlation
plot(1000*time_cor, squeeze(meas_cor(:,nIR,band)) ,'-', 'color', cMap2(:, 1), 'LineWidth',3)
% expected correlation from the noise
plot(1000*time_cor, squeeze(r_snr(:, nIR, band)), 'color', cMap2(:, 2), 'LineWidth',2)
%model prediction
plot( 1000*time_cor, pred_cor_resample(:,nIR,  band) ,'-.',  'color', cMap2(:, 3), 'LineWidth',2)
% plot exptected + prediction
plot(  1000*time_cor, squeeze(r_snr(:, nIR, band)).*pred_cor_resample(:,nIR, band),'--' ,'color', cMap2(:, 4), 'LineWidth',2)



ylim([0.9 1.01])
xlim(1000*[0 0.55]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation', 'Interpreter','latex', 'FontSize',12)
set(gca, 'ytick', .9:.02:1, 'xtick', 0:100:600, 'FontSize',12)
lgd  = legend('$\rho_{x_i, x_j}$', '$\widehat{\rho}_{x_i, x_j}$', '$\rho_{1r}$', '$\widehat{\rho}_{x_i, x_j} \times \rho_{1r}$', 'location', 'south', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
box on
grid on
f.Position(end) = 270;
%% print figure 3
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_different_correlations','-dpdf','-r0')

%% interpolate the correlation 
for b = 1:bandIt
    for i =1: numRIR-1
        corr_temp(:, i, b) = interp1( pred_toastd ./ (sqrt(2)*divergence(i,b)), pred_cor(:, b) , time_cor);
    end
end
%% plot Arni - different bands
offset = 0;
f = figure(1); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band = [10 19];
for b = band
    for i = 1
        plot( 1000*time_cor, squeeze(meas_cor(:,1)) + offset(end),'-', 'color', cMap2(:, ceil(b/10)), 'LineWidth',3)
        
%         plot(  1000*time_cor, median(exp_corr(:, [refIR, i+1], b),2).*corr_temp(:, i, b),'-' ,'color', cMap2(:, ceil(b/10)+2), 'LineWidth',1.5)

    end
end

ylim([0 1.01])
xlim(1000*[0 0.65]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation', 'Interpreter','latex', 'FontSize',12)
lgd  = legend('$\rho_{x_i, x_j}$, 10 kHz', '$\widehat{\rho}_{x_i, x_j} \times \rho_{1r}$, 10 kHz', '$\rho_{x_i, x_j}$, 19 kHz', '$\widehat{\rho}_{x_i, x_j} \times \rho_{1r}$, 19 kHz','location', 'southwest', 'interpreter', 'latex', 'fontsize', 11, 'numcolumns', 1);


box on
grid on
f.Position(end) = 320;

%% print figure
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_correlation_Jan2024','-dpdf','-r0')

%% get the correlation coefficient between the calculated correlation curves and the fitted curves
clear new_cor

cor_PCC = zeros(2,2,4,bandIt);
for i = 1:numRIR-1
    for it = 1: bandIt
        firstI = find(mask(:, it) == 1, 1, 'first');
        firstI = firstI+500;
        corL = find(mask(firstI+1:end, it) == 0, 1, 'first');
        nc_len = corL;
       
        cor_PCC (:, :, i, it) = corrcoef([meas_cor(firstI:firstI+nc_len-1,i, it), corr_temp(firstI:firstI+nc_len-1, i, it)]);
    end
end
%% plot the correlation coefficient for the fits
f = figure(7); clf; hold on

for i =1:numRIR-1 %[1 6 10]
    plot(bandCenters, squeeze(cor_PCC(1,2,i, :)), '-.o', 'color', cMap2(:, i), 'LineWidth',1, 'MarkerFaceColor',cMap2(:, i))
end 
% yline(0.7, 'k--', 'HandleVisibility','off')
ylim([0.75 1])
 xlim(1000*[0.5 19.5]);
xlabel('Center frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel([{'Model-signal correlation'}], 'Interpreter','latex', 'FontSize',12)
set(gca,  'xtick', 1000*[1:2:20], 'xticklabel', 1:2:20,  'FontSize',12)
% set(gca, 'ytick',[0.75:0.1:1, 1])
% set(gca, 'ytick',[-0.5 0 0.5 0.7 1])
box on
f.Position(end) = 250;
% 
% lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'southeast', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);

lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'south', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 4);
% lgd  = legend('5 min', '30 min', '50 min',  'location', 'northwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 3);

lgd.Title.String = [ {'Time between measurements'}];%$j-i$';
%% print the correlation coefficient figure
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'corr_coef_Arni','-dpdf','-r0')
% print(f,'corr_coef_Arni_FA','-dpdf','-r0')


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
        F = fb(numBands);

        % l1 loss
        loss_fun = @(volatility) sum(abs(correlationModel(F,T,exp(volatility)).*snr - cor));

        % do the fitting in log(volatility) for better scaling
        volatilityMin = -50;
        volatilityMax = -3;
        options = optimset('TolX',0.01);
        vol = exp(fminbnd(loss_fun,volatilityMin,volatilityMax,options));

        volatility(itIR,itBands) = vol;

        figure; hold on;
        plot(T,cor)
        plot(T,correlationModel(F,T,vol).*snr);
        
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