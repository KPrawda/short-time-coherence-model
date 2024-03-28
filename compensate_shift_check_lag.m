%% compensate for temperature drift with Postma-Katz method
clear all; close all; clc
set(groot,'defaultAxesTickLabelInterpreter','latex');
%%
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Mesurements Arni\IR_Arni_upload\IR_Arni_upload';
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR';
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\RIRs'
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Forum Acusticum 2023\RIRs\'
%%
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

folder='C:\Users\prawdak1\Dropbox (Aalto)\Projects\Forum Acusticum 2023\RIRs\';
audio_files=dir(fullfile(folder,'*.wav'));

file = cell(numel(audio_files),1);
for k=1:numel(audio_files)
  filename=audio_files(k).name;

      file{k} = filename;     

end
%% set the first RIR as a reference RIR
[ir(:, 1), fs] =  audioread(file{1});


%% upsampling and windowing parameters
upFactor = 10; %upsample by this factor
winLen = 0.01*fs*upFactor; % window length of 10 ms 
wL = 0.01*fs; % not up-sampled window length
%% set the first RIR as a reference RIR
[ir(:, 1), fs] =  audioread(file{1});
ir_upsampled(:, 1) = interp(ir(:, 1), upFactor);
referenceIR=1; % set RIR as the reference 

newPoints = (0:upFactor*length(ir)-1)./(upFactor*fs);


%% read the RIRs
folder_in='C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR\Up-to-date code\Resampled_RIRs\';
% clear ir
% 1-1293 rec1 and rec 2, the rest rec 3
it = 2;
for ind = 1137:1293 %% read the IRs: IR_numClosed_50_numComb_5000_mic_1_sweep_x.wav    21:24%
    tempIR =  audioread(file{ind});
    ir(:, it) = tempIR(1:length(ir));
    if sum(ir(:, it)) ==0
        continue
    end

   [value, location] = max((ir));
    if diff(location) > 0
        ir(:, 2) = [ir(diff(location)+1:end, 2); ir(end -diff(location)+1:end, 2)];
    elseif diff(location) < 0
        ir(:, 1) = [ir(diff(location)+1:end, 1); ir(end -diff(location)+1:end, 1)];
    end




%% upsample RIRs
%     ir_upsampled = zeros(length(ir)*upFactor,2);
    ir_upsampled(:, it) = interp(ir(:, it), upFactor);
   

%% evaluate the lag between the reference RIR and the time-stretched RIR
 
    [maxLag, snrDB] = getCorrLag(ir_upsampled(:, [referenceIR, it]),winLen, referenceIR); % get the correlation lag
%     [~, windMin] =  min(abs(maxLag(:, end)));%find(maxLag(:, end) == 0, 1, "first");
    [~, loc] = min(abs(snrDB-15)); % restrict the analysis (slope fitting) only to the range with high enough SNR
    
    winInd = round(loc(1)./winLen); % index of the last useful window 

%% evaluate the slope of the max lag for time-stretching factor
    xx = (1:max(winInd))*wL/fs; % take the windows with high enough snr
    xx = xx.';
    y = maxLag(1:max(winInd), end); % take the lag values for which snr is high enough
    xx = xx.*ones(1,size(y, 2));  % make as many vectors as we have RIRs
    fun = @(m)(sum(abs(xx*m-y))); % linear function for estimating the slope
    m(it) = fminsearch(fun, 0);                % slope estimation
    

%% resample
 m
     if m(it)~= 0
%          temp_ = resample(rir,upFactor*fs, round(upFactor*fs-m(it)));
         ir_res= resampleGiantFFT(ir_upsampled(:, it), round(upFactor*fs-m(it)),upFactor*fs);
         
         ir_down = downsample(ir_res, upFactor);
         
         
     else 
         ir_down = zeros(size(ir,1), 1);
         ir_down(1:length(ir)) = ir(:, it);
     end
     fil = file{ind};
     filename = strcat(fil(1:end-4), '_resampled.wav');
     full_filename = [folder_in, filename];
     audiowrite(full_filename, ir_down, fs, 'BitsPerSample',32)
end
%%
% temp_3 = resampleGiantFFT(rir, upFactor*fs,fs);

%%
len = min(length(ir_upsampled), length(ir_res));

 [maxLag_, snrDB_] = getCorrLag([ir_upsampled(1:len, referenceIR), ir_res(1:len)],winLen, referenceIR);%, fs, upFactor);


%% plot the RIRs
oldPoints = (0:length(ir)-1)./fs;
f = figure(2); clf

subplot(2,1,1); hold on; box on; grid on
plot(1000*oldPoints, db(ir(:, 1)./max(abs(ir(:, 1)))), 'k', 'LineWidth',2)
plot(1000*oldPoints, db(ir(:, end)/max(abs(ir(:, end)))), 'color',col1,  'LineWidth',1)
ylim([-100 1])
xlim([0 1200])
%% plot the time lag vs the RIR difference in dB



f = figure(1); clf; hold on

% plot the lag
subplot(2,1,1); hold on
plot(1000*[0:winInd]*winLen./(upFactor*fs), maxLag(1:winInd+1, 2), 'Color', cMap2(:, 1), 'LineWidth',1.2) % between ref RIR(1) and RIR 2
plot(1000*[0:winInd]*winLen./(upFactor*fs), maxLag(1:winInd+1, 4), 'Color', cMap2(:, 3), 'LineWidth',1.2) % between ref RIR(1) and RIR 4

xlim([0 1200])
ylim([- 3.1 2.5])

set(gca,  'FontSize',12)
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Time lag (samples)', 'Interpreter','latex', 'FontSize',12)
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

%% sliding correlation
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

%% resampling using giant FFTs
function outSignal = resampleGiantFFT(inSignal, oldFs, newFs)

resampling_factor = oldFs/newFs;
if resampling_factor == 1
    outSignal = inSignal;
else 
    sigLen = ceil(max(oldFs/newFs, newFs/oldFs)*length(inSignal));
    temp_signal = zeros(sigLen, 1);
    temp_signal(1:length(inSignal)) = inSignal;
    if resampling_factor < 1
        sig_FFT = fft(temp_signal, oldFs);
        res_FFT = [sig_FFT(1: oldFs/2); ones(newFs-oldFs,1)*10^-16; sig_FFT(end - oldFs/2 +1:end)]; 
        outSignal = real(ifft(res_FFT));
    elseif resampling_factor > 1
        sig_FFT = fft(temp_signal, oldFs);
        res_FFT = [sig_FFT(1: newFs/2); sig_FFT(end - newFs/2 +1:end)]; 
        outSignal = real(ifft(res_FFT));
    end


end

end
