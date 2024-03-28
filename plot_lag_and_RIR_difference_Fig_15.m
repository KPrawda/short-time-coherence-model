%% calculate the time lag and RIR difference and compute the compensation factor for temperature drift with Postma-Katz method
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024
%% housekeeping 
clear all; close all; clc
set(groot,'defaultAxesTickLabelInterpreter','latex');
addpath('.\RIRs')
%% colors
numPlots =4 ;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cMap3 = [linspace(col1(1),1,numPlots); linspace(col1(2),1,numPlots); linspace(col1(3),1,numPlots)];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;
%% read the corresponding IR from the Arni database
folder='.\RIRs';
audio_files=dir(fullfile(folder,'*.wav'));

file = cell(numel(audio_files),1);
for k=1:numel(audio_files)
    filename=audio_files(k).name;
    file{k} = filename;     
end
%% read the RIRs
clear ir
li = 1;
for it = 6:10 %% read the IRs: IR_numClosed_50_numComb_5000_mic_1_sweep_x.wav    21:24%
    [ir(:, li), fs] =  audioread(file{it});
    li = li+1;
end
%% delete the pre-delay
[~, loc_max] = max(abs(ir)); 
if loc_max>1
    for it = 1:li-1
        ir(:, it) = [ir( loc_max(it):end, it); zeros(loc_max(it)-1, 1)];   
    end
end
%% upsampling and windowing parameters
upFactor = 10; %upsample by this factor
winLen = 0.01*fs*upFactor; % window length of 10 ms 
wL = 0.01*fs; % not up-sampled window length
%% upsample RIRs
numRIR = min(size(ir)); % number of RIRs
lenRIR = max(size(ir)); % RIR lengths (samples)

ir_upsampled = zeros(lenRIR*upFactor, numRIR); % initialize the upsampled RIR

for it = 1:numRIR
    ir_upsampled(:, it) = interp(ir(:, it), upFactor);
end
%% evaluate the lag between the reference RIR and the time-stretched RIR
referenceIR=1; % which RIR is the reference 
[maxLag, snrDB] = getCorrLag(ir_upsampled,winLen, referenceIR); % get the correlation lag
[~, loc] = min(abs(snrDB-5)); % restrict the analysis (slope fitting) only to the range with high enough SNR
winInd = round(loc./winLen); % index of the last useful window 
%% evaluate the slope of the max lag for time-stretching factor
xx = (1:max(winInd))*wL/fs; % take the windows with high enough snr
xx = xx.';
y = maxLag(1:max(winInd), :); % take the lag values for which snr is high enough
xx = xx.*ones(1,size(y, 2));  % make as many vectors as we have RIRs
m = zeros(1, size(y, 2));     % initialize slope values
for i = 1:size(y, 2)
    fun = @(m)(sum(abs(xx(:, i)*m-y(:, i)))); % linear function for estimating the slope
    m(i) = fminsearch(fun, 0);                % slope estimation
end
%% plot the time lag vs the RIR difference in dB, figure 15

f = figure(1); clf; hold on

% plot the lag
subplot(2,1,1); hold on
plot(1000*[0:winInd]*winLen./(upFactor*fs), (10^6)*maxLag(1:winInd+1, 2)./(upFactor*fs), 'Color', cMap2(:, 1), 'LineWidth',1.2) % between ref RIR(1) and RIR 2
plot(1000*[0:winInd]*winLen./(upFactor*fs), (10^6)*maxLag(1:winInd+1, 4)./(upFactor*fs), 'Color', cMap2(:, 3), 'LineWidth',1.2) % between ref RIR(1) and RIR 4

xlim([0 1200])
% ylim([- 3.1 2.5])

set(gca,  'FontSize',12)
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Time lag ($\mu$s)', 'Interpreter','latex', 'FontSize',12)
box on


lgd  = legend('5 s',  '15 s',  'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
lgd.Title.String = [ {'Time between measurements'}];%$j-i$';

% plot thr RIR difference in dB
time = [0:length(ir)-1]./fs; % time vector for the original RIRs

subplot(2,1,2); hold on
plot(1000*time, db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 4)./max(abs(ir(:, 4)))), 'Color', cMap2(:, 3), 'LineWidth',1 )
plot(1000*time, db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 2)./max(abs(ir(:, 2)))), 'Color', cMap2(:, 1), 'LineWidth',1 )

xlim([0 1200])
ylim([-80 -35])
box on
set(gca,  'FontSize',12)
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('RIR difference (dB)', 'Interpreter','latex', 'FontSize',12)
lgd  = legend('5 s',  '15 s',  'location', 'northeast', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
lgd.Title.String = [ {'Time between measurements'}];%$j-i$';
%% print the figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'IR_difference','-dpdf','-r0')
%% estimate the differences between RIR in dB
% noise floor level of the difference of RIRs
median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end,2)./max(abs(ir(:, 2)))))
median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end, 4)./max(abs(ir(:, 4)))))

% max level of the difference of rirs
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 2)./max(abs(ir(:, 2)))))
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 4)./max(abs(ir(:, 4)))))

% difference between the two
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 2)./max(abs(ir(:, 2))))) - median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end,2)./max(abs(ir(:, 2)))))
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 4)./max(abs(ir(:, 4))))) - median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end, 4)./max(abs(ir(:, 4)))))
%%  calculate lag
function [maxLag, snrDB] = getCorrLag(ir,winLen, referenceIR)

num_rirs = min(size(ir)); % numbers of RIRs (channels)

L = floor(size(ir,1)/winLen); % number of window to evaluate the lag
maxCor = zeros(L, num_rirs);
maxLag = zeros(L, num_rirs);

% calculate the lag for each window
for it = 1:num_rirs % for each channel    
    for n = 1:L % for each window
        [c, lags] = xcorr(ir((winLen*(n-1)+1: winLen*n), [ referenceIR, it]), 'normalized'); % calculate cross-correlation between the ref IR and IR in question
        [maxCor(n, it), tempLag] = max(abs(c(:, 2))); % find the max value of the cross-correlation
        maxLag(n, it) = lags(tempLag); % find the lag of max cross-corr
    end
end
win = rectwin(winLen);%hann(winLen) is much smoother than rectangular window, but we use rectangular anyway
win = win / sum(win); % normalize the window 

for it = 1:num_rirs % for each channel
    energy(:,it) = conv(ir(:,it).^2,win,'valid'); % calculate energy for the whole signal
end
energy(energy(:, 1)==0, :)=[];
enoise =  median(energy(end-2*winLen:end-winLen , :)); % estimate the noise energy by taking the energy of the last window
SNR = abs(energy-enoise)./enoise; % signal-to-noise ratio
snrDB = 10*log10(SNR);

end
