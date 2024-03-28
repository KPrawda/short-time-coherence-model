%% compensate for temperature drift with Postma-Katz method
clear all; close all; clc
set(groot,'defaultAxesTickLabelInterpreter','latex');
%%
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Mesurements Arni\IR_Arni_upload\IR_Arni_upload';
addpath 'C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR';
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

folder='C:\Users\prawdak1\Dropbox (Aalto)\Projects\Time-variance\TimeVaryingCorrelationRIR';
audio_files=dir(fullfile(folder,'*.wav'));

file = cell(numel(audio_files),1);
for k=1:numel(audio_files)
  filename=audio_files(k).name;

      file{k} = filename;     

end
file = file';
fs = 44100;
%%
clear ir
li = 1;
for it = 14:18 %% read the two IRs
    
    ir(:, li) =  audioread(file{it});
    li = li+1;
end

[~, loc_max] = max(abs(ir));
for it = 1:li-1
    ir(:, it) = [ir( loc_max(it):end, it); zeros(loc_max(it)-1, 1)];

end
ir = ir(1:end-max(loc_max), :);
%% upsampling and windowing parameters
upFactor = 10;
winLen = 0.01*fs*upFactor; % window length of 100 ms 
wL = 0.01*fs;

%% evaluate the lag between the reference RIR and the time-stretched RIR
referenceIR=1; % which RIR is the reference 
[maxLag, snrDB] = getCorrLag(ir,winLen, referenceIR, fs, upFactor); % get the correlation lag
[~, loc] = min(abs(snrDB-5)); % restrict the analysis (slope fitting) only to the range with high enough SNR
loc = loc*upFactor;
winInd = round(loc./winLen); % index of the last useful window 

%% evaluate the slope of the max lag for time-stretching factor
xx = (1:max(winInd))*wL/fs;
xx = xx.';
y = maxLag(1:max(winInd), :);
xx = xx.*ones(1,size(y, 2));
m = zeros(1, size(y, 2));
for i = 1:size(y, 2)
    fun = @(m)(sum(abs(xx(:, i)*m-y(:, i))));
    m(i) = fminsearch(fun, 0);
end
%%
oldPoints = (0:length(ir)-1)./fs;
f = figure(2); clf

subplot(2,1,1); hold on; box on; grid on
plot(1000*oldPoints, ir(:, 1), 'k', 'LineWidth',2)
plot(1000*oldPoints, ir(:, end), 'color',col1,  'LineWidth',1)
ylim([-1 1])

%% plot the time lag - ARni
time = [0:length(ir)-1]./fs;
f = figure(1); clf; hold on
subplot(2,1,1); hold on
plot(1000*[0:winInd]*winLen./(upFactor*fs), maxLag(1:winInd+1, 2), 'Color', cMap2(:, 1), 'LineWidth',1.2)
plot(1000*[0:winInd]*winLen./(upFactor*fs), maxLag(1:winInd+1, 4), 'Color', cMap2(:, 3), 'LineWidth',1.2)

xlim([0 1200])
ylim([- 3.1 2.5])
box on
% grid on
set(gca,  'FontSize',12)
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Time lag (samples)', 'Interpreter','latex', 'FontSize',12)
% f.Position(end) = 200;

lgd  = legend('5 s',  '15 s',  'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
lgd.Title.String = [ {'Time between measurements'}];%$j-i$';

subplot(2,1,2); hold on

plot(1000*time, db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 4)./max(abs(ir(:, 4)))), 'Color', cMap2(:, 3), 'LineWidth',1 )
plot(1000*time, db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 2)./max(abs(ir(:, 2)))), 'Color', cMap2(:, 1), 'LineWidth',1 )

xlim([0 1200])
ylim([-80 -35])
box on
% grid on
set(gca,  'FontSize',12)
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('RIR difference (dB)', 'Interpreter','latex', 'FontSize',12)

%%
% noise floor level of the difference of RIRs
median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end,2)./max(abs(ir(:, 2)))))
median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end, 4)./max(abs(ir(:, 4)))))

% max level of the difference of rirs
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 2)./max(abs(ir(:, 2)))))
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 4)./max(abs(ir(:, 4)))))

% difference between the two
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 2)./max(abs(ir(:, 2))))) - median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end,2)./max(abs(ir(:, 2)))))
max(db(ir(:, 1)./max(abs(ir(:, 1)))-ir(:, 4)./max(abs(ir(:, 4))))) - median(db(ir(fs:end, 1)./max(abs(ir(:, 1)))-ir(fs:end, 4)./max(abs(ir(:, 4)))))
%% resample

newPoints = (0:upFactor*length(ir)-1)./(upFactor*fs);
for it = 1:size(y, 2)
 rir = interp(ir(:, it), upFactor);
 round(upFactor*fs-m(it))
     if m(it)~= 0
         temp_ = resample(rir, upFactor*fs, round(upFactor*fs-m(it)));
         ir_res(1:length(temp_), it) = temp_;
         ir_down(1:length(temp_)/upFactor+1, it) = downsample(ir_res(:, it), upFactor);
     else 
         ir_down(:, it) = zeros(size(ir,1), 1);
         ir_down(1:length(ir), it) = ir(:, it);
     end
end

%%
[maxLag_, snrDB_] = getCorrLag([ir(:, referenceIR), ir_down(1:length(ir), referenceIR+1:end)],winLen, referenceIR, fs, upFactor);



%% sliding correlation
function [maxLag, snrDB] = getCorrLag(ir,winLen, referenceIR, fs, upFactor)
% referenceIR = 1;
num_rirs = size(ir,2)
% winLen = 0.01*fs*upFactor;
% oldPoints = (0:length(ir)-1)./fs;
% newPoints = (0:upFactor*length(ir)-1)./(upFactor*fs);
% rir = interp1(oldPoints, ir, newPoints);
for it = 1:num_rirs % for each channel
rir(:, it) = interp(ir(:, it), upFactor);
end

L = floor(size(rir,1)/winLen);
maxCor = zeros(L, num_rirs);
maxLag = zeros(L, num_rirs);
for it = 1:num_rirs % for each channel
    
    for n = 1:L % for each window
        [c, lags] = xcorr(rir((winLen*(n-1)+1: winLen*n), [ referenceIR, it]), 'normalized');
        [maxCor(n, it), tempLag] = max(abs(c(:, 2)));
        maxLag(n, it) = lags(tempLag);
    end
end
win = hann(winLen); % much smoother than rectangular window
win = win / sum(win);

for it = 1:num_rirs % for each channel
    energy(:,it) = conv(ir(:,it).^2,win,'valid');
end
enoise =  energy(end-winLen, :);
SNR = abs(energy-enoise)./enoise;
snrDB = 10*log10(SNR);

end