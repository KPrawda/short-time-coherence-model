% Sebastian J. Schlecht, Sunday, 04 December 2022
clear; clc; close all;
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\RIRs'
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Forum Acusticum 2023\RIRs'
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\Resampled_RIRs\'
%%
% in Arni for mic 1 to 5
% direct_delay = [560   625   334   224   187]; % samples
switch 1
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
        numRIR = 4;
        filename = 'flutter4_setting2_IR_channel_4_sweep_%d.wav';
        direct_delay = 560;
        for it = 2:numRIR+1
            [rir(:,it-1),fs] = audioread(sprintf(filename,it));
        end
        [~,onset] = max(rir(:,2));
        fit_onset = floor(onset) - direct_delay;
        rir = rir(fit_onset+500:end, :);
        minDB = 30;
    case 4
        numRIR = 2;
        filename = 'IR_Synthetic_%d.wav';
        direct_delay = 0;
        fit_onset = 1;
        minDB = 0;
    case 5
        numRIR = 55;
        filename = 'rec3_IR_ch_8_%d.wav';
        direct_delay = 0;
        fit_onset = 1;
        minDB = 30;
        its = 1:55;
        for it = its%[1,53: 55]%1:numRIR
            [rir(:,it),fs] = audioread(sprintf(filename,it));
        
        end
        rir = rir(:, its);
        for it = 1:numRIR
            if rir(1, it)< -0.003
                rir(:, it) = [rir(3:end,it);rir(end,it);rir(end,it)];
            end
        end
    case 6    
        filename = 'rec1_IR_ch_8_';
        minDB = 30;
        its = [ 60:60:600];
        numRIR = length(its)+1;
        [rir(:, 1), fs] = audioread('rec1_IR_ch_8_1.wav');
        for it = its
            temp = audioread(sprintf('%s%d%s',filename,it, '_resampled.wav'));
            rir(1:length(temp),it) = temp;
        
        end
        rir = rir(:, [1,its]);
end


% rir = rir(fit_onset:end,:); % truncate to time of sound emittance
% rir = circshift(rir,-1,2); % reorder the sweeps​
% ​
%% bandpass filtering and correlation estimation
winLen = 2^10;
bandCenters =  (1:19)*1000; % Hz
refIR = 1;
for bandIt = 1:numel(bandCenters)
    band_freq = [-500 500]+bandCenters(bandIt);
    rir_band = bandpass(rir,band_freq,fs);

    % compute correlation
    [cor, energy] = slidingCorrelation(rir_band, winLen);

    meas_cor(:,:,bandIt) = cor(:,2:end);

    % find valid part
    energyDB = db(mean(energy,2));
    mask(:,bandIt) = energyDB > min(energyDB) + minDB ;
    mask(2*fs:end,bandIt) = 0;
    %     mask(1:(winLen/2 + direct_delay),:) = 0;

    % expected correlation based on SNR
%     eNoise= (energy(round(0.7*length(energy)):round(0.8*length(energy)),: ));
    eNoise= (energy(round(0.8*fs):round(0.9*fs),: ));
    eN = median(eNoise, 1);
    snr = (energy-eN)./eN;
    
    snrdB = 10*log10(abs(snr));
%     [~, locs] = min(abs(snrdB(1: 20000, :)));
    
    exp_corr(:,:,bandIt)  = snr./(1+snr);
end
time_cor = (1:size(cor,1)).'/fs; % seconds

%% towards analytical model
[pred_cor,pred_toastd] = singlePulseCorrelation(bandCenters,fs);

%% find divergence
% load("Arni_div_FA_new.mat")
[divergence] = findDivergence(pred_toastd, pred_cor, time_cor, meas_cor, mask);
% divergence = v_s(:, its);
% divergence = divergence';
% findDivergence(pred_toastd, pred_cor, time_cor, meas_cor(:, :, end), mask(:,  end));
% save('div_sim.mat', 'divergence')
%% plot
figure; hold on; grid on;
plot(bandCenters,divergence')
xlabel('Frequency Band (Hz)');
ylabel('Divergence (seconds)');

% plots only first two RIR correlation
figure; hold on; grid on;
offset = 0*(1:size(pred_cor,3));
plot( time_cor, squeeze(meas_cor(:,3,:)) + offset)
set(gca,'ColorOrderIndex',1);
plot( pred_toastd ./ (sqrt(2)*divergence(3,:)), pred_cor + offset)
xlim([0 max(time_cor)])


%% colors
numPlots = numRIR-1;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;
%% plot Arni - measured, modeled, expected and so on 
offset = 0;
f=figure(3); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band = 10; % which band to plot
nIR = 5; % which RIR to plot
%measured correlation
plot(1000*time_cor, squeeze(meas_cor(:,nIR,band)) + offset(end),'-', 'color', cMap2(:, 1), 'LineWidth',3)
% expected correlation from the noise
plot(1000*time_cor, median(exp_corr(:, [refIR, nIR], band),2), 'color', cMap2(:, 2), 'LineWidth',2)
%model prediction
plot( 1000*pred_toastd ./ (sqrt(2)*divergence(nIR,band)), pred_cor(:, band) + offset(end),'-.',  'color', cMap2(:, 3), 'LineWidth',2)
% extrapolate the model prediction
corr_temp = interp1( pred_toastd ./ (sqrt(2)*divergence(nIR,band)), pred_cor(:, band) + offset(end), time_cor);
% plot exptected + prediction
plot(  1000*time_cor, median(exp_corr(:, [refIR, nIR+1], band),2).*corr_temp,'--' ,'color', cMap2(:, 4), 'LineWidth',2)



ylim([0.9 1.01])
xlim(1000*[0 0.55]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation', 'Interpreter','latex', 'FontSize',12)
set(gca, 'ytick', .9:.02:1, 'xtick', 0:100:600, 'FontSize',12)
% % lgd  = legend('$1$', '$2$', '$3$', '$4$', 'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
% % lgd.Title.String = '$j-i$';
lgd  = legend('$\rho_{x_i, x_j}$', '$\widehat{\rho}_{x_i, x_j}$', '$\rho_{1r}$', '$\widehat{\rho}_{x_i, x_j} \times \rho_{1r}$', 'location', 'south', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
% % lgd.Title.String = [ {'Time between'},{' measurements'}];%$j-i$';
box on
grid on
f.Position(end) = 270;
% plot(1000*time_cor, mask(:, end), 'k:')
%% print figure 3
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
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
    for i = 2
        plot( 1000*time_cor, squeeze(meas_cor(:,i,b)) + offset(end),'-', 'color', cMap2(:, ceil(b/10)), 'LineWidth',3)
        
        plot(  1000*time_cor, median(exp_corr(:, [refIR, i+1], b),2).*corr_temp(:, i, b),'-' ,'color', cMap2(:, ceil(b/10)+2), 'LineWidth',1.5)

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
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Arni_correlation_Jan2024','-dpdf','-r0')
%% plot Arni - different time separations
offset = 0;
f = figure(1); clf; hold on

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
band = [18];
for b = band
    for i = [1 6 10]
        plot( 1000*time_cor, squeeze(meas_cor(:,i,b)) + offset(end),'-', 'color', cMap2(:, i), 'LineWidth',2)
        
        plot(  1000*time_cor, median(exp_corr(:, [refIR, i+1], b),2).*corr_temp(:, i, b),'--' ,'color', cMap2(:, i), 'LineWidth',1.5, 'HandleVisibility','off')

    end
end

ylim([0.7 1.01])
xlim(1000*[0 0.2]);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation', 'Interpreter','latex', 'FontSize',12)
% set(gca, 'ytick', .0:.2:1, 'xtick', 0:100:700, 'FontSize',12)

lgd  = legend('5 min', '30 min', '50 min',  'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 3);

lgd.Title.String = [ {'Time between measurements'}];%$j-i$';
box on
grid on
f.Position(end) = 320;
%% print figure 3
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Arni_FA_correlation','-dpdf','-r0')


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
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'corr_coef_Arni','-dpdf','-r0')
% print(f,'corr_coef_Arni_FA','-dpdf','-r0')
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
win = rectwin(winLen); % hann much smoother than rectangular window
win = win / sum(win);

for it = 1:num_rirs % for each channel
    xCov(:,it) = conv(ir(:,it) .* ir(:,referenceIR),win,'valid');
    energy(:,it) = conv(ir(:,it).^2,win,'valid');
    diffEnergy(:,it) = conv( (ir(:,it) - ir(:,referenceIR)).^2,win,'valid');
end
xCor = xCov ./ sqrt(energy) ./ sqrt(energy(:,referenceIR));

end