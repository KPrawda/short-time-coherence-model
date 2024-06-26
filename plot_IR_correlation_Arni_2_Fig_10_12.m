%% plot figures 10 and 12 from the manuscript
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024

% Sebastian J. Schlecht, Sunday, 04 December 2022
% new version 2024-01-26
% updated by K. Prawda 12.02.2024 and 7.06.2024
% PLOT ONLY THE CORRELATION FOR ARNI 2 DATASET
%% housekeeping
clear; clc; close all;
addpath(fullfile(pwd, 'Resampled_RIRs'));
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%% Load measurements
referenceRIR = 1;

filename = 'rec1_IR_ch_8_';
minDB = 30;
its =  60:60:600;
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

    [cor, energy, r_snr_temp, e_ref, e_rir] = slidingCorrelation(rir_band_ref, rir_band_other, winLen);
    
    r_snr(:,:,bandIt) = r_snr_temp;
    meas_energy(:,:,bandIt) = energy;
    meas_cor(:,:,bandIt) = cor;
    meas_energy_ref(:,:,bandIt) = e_ref;

    % only estimate volatility for the part of IR with sufficiently
    % high energy
    energyDB = db(mean(energy,2));
    mask(:, bandIt) = energyDB > min(energyDB) + minDB ;
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

%% plot Arni 2 - different time separation, measured and modeled coherence, figure 10
f = figure(2); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 

b =18; % which band to plot

 for i = [1 6 10] % which RIRs to plot
    plot( 1000* time_cor, squeeze(meas_cor(:,i,b)).^2 ,'-', 'color', cMap2(:,i), 'LineWidth',2)
    plot(  1000*time_cor, (squeeze(r_snr(:, i, b)).*pred_cor_resample(:,i, b)).^2,'--' ,'color', cMap2(:, i), 'LineWidth',2, 'HandleVisibility','off')   
 end    


ylim([0.7 1.03])
xlim([0 200])
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Coherence', 'Interpreter','latex', 'FontSize',12)
set(gca, 'XTick', 0:50:200, 'fontsize', 12)
lgd  = legend('5 min', '30 min', '50 min',  'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 3);
lgd.Title.String = [ {'Time between measurements'}];%$j-i$';
box on
grid on
f.Position(end) = 320;

%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'Arni_FA_correlation','-dpdf','-r0')

%% get the correlation coefficient between the calculated correlation curves and the fitted curves
cor_PCC = zeros(2,2,numRIR-1,numel(bandCenters));
for i = 1:numRIR-1
    for it = 1: numel(bandCenters)
        firstI = find(mask(:, it) == 1, 1, 'first');
        corL = find(mask(firstI+1:end, it) == 0, 1, 'first');
        nc_len = corL;
        corr_temp =  (squeeze(r_snr(:, i, it)).*pred_cor_resample(:,i, it)).^2;
       
        cor_PCC (:, :, i, it) = corrcoef([meas_cor(firstI:firstI+nc_len-1,i, it).^2, corr_temp(firstI:firstI+nc_len-1)]);
    end
end
%% plot the model-signal correlation coefficient for the fits, figure 12
f = figure(3); clf; hold on

for i =[1 6 10]
    plot(bandCenters, squeeze(cor_PCC(1,2,i, :)), '-.o', 'color', cMap2(:, i), 'LineWidth',1, 'MarkerFaceColor',cMap2(:, i))
end 

ylim([0.89 1.01])
xlim(1000*[0.5 19.5]);
xlabel('Center frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('Model-signal correlation', 'Interpreter','latex', 'FontSize',12)
set(gca,  'xtick', 1000*[1:2:20], 'xticklabel', 1:2:20,  'FontSize',12)
box on
f.Position(end) = 250;
lgd  = legend('5 min', '30 min', '50 min',  'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 3);
lgd.Title.String = [ {'Time between measurements'}];
%% print the correlation coefficient figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'corr_coef_Arni_FA','-dpdf','-r0')
