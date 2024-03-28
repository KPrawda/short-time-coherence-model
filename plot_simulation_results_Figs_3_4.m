%% plot short-time coherence simulation results - Figures 3 and 4 from the manuscript
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024

% Sebastian J. Schlecht, Tuesday, 16. January 2024
% test time-frequency correlation model based on normal TDoA distribution
%% housekeeping
clear; clc; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%%
fs = 48000*1;
lenSeconds = 2;
lenSamples = fs*lenSeconds;

% generate echo times
num_echoes = 2*10^4;
mean_echo_times = rand(num_echoes,1) * lenSeconds;

% generate variation
volatility = 10^-6*5; % TOASTD increase per second.
num_rirs = 2;

toaVariation = sqrt(mean_echo_times) .* volatility .* randn(num_echoes,num_rirs); % Gaussian random walk
toaVariation(:,1) = 0; % no variation on the reference RIR
var_echo_times = mean_echo_times + toaVariation;

var_echo_times = clip(var_echo_times, [0 lenSeconds]);

% generate signs
echo_amps = 2*round(rand(num_echoes,1))-1;

% synthesis with fractional delays and amplitudes
for it = 1:num_rirs
    rir(:,it) = synthRIR(var_echo_times(:,it), echo_amps, fs, lenSeconds);
end

% audiowrite
rirPeakValue = max(abs(rir(:)));
for it = 1:size(rir,2)
    filename = 'IR_Synthetic_%d.wav';
    audiowrite(sprintf(filename,it),rir(:,it) ./ rirPeakValue * 0.99,fs);
end

%% time-frequency correlation based on STFT
clear s
for it = 1:size(rir,2)
    [s(:,:,it),F,T] = stft(rir(:,it),fs,Window=kaiser(256,5),OverlapLength=240,FFTLength=1024);
end

s = s(F >= 0,:,:); % remove negative frequencies
F = F(F >= 0);

% compute correlation 
freqWin = 40; 
timeWin = 180;

Energy = abs(s).^2;
Energy = movsum(Energy,freqWin,1);
Energy = movsum(Energy,timeWin,2);
Energy = mean(Energy,3); % average channels

Cov = s(:,:,1) .* conj(s(:,:,2));
Cov = movsum(Cov,freqWin,1);
Cov = movsum(Cov,timeWin,2);

Cor = abs(Cov) ./ Energy;

figure;
surf(T,F,Cor,'EdgeColor','none')
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
view([0 90])
colorbar



%% colors
numPlots = 20;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];
cMap2 = cMap2';
col2 = [113, 62, 90]./255;
%% the model of time-frequency correlation is based on the Gaussian TDoA
magicFactor = 20; % might be influenced by the STFT parameters
CovModel = exp( - magicFactor * (F .* sqrt(T).'*volatility).^2 );
%% plot figure 3 
% Freqs = (1:20)'*1000;
T_ = [0.1:0.05:1].';
CovMod = exp( - magicFactor * (F .* sqrt(T_).'*volatility).^2 );

f = figure(5);clf; hold on; grid on; box on
plot(F./1000,CovMod.^2);
% xlim([0 fs/2]./1000)
xlim([0 20])
colororder(cMap2(1:length(T_), :))

xlabel('Frequency (kHz)', 'Interpreter','latex')
ylabel('$R^2_\vartheta (\tau, f)$', 'Interpreter','latex')
set(gca, 'FontSize', 12, 'XTick', [0:5:20, 24])
clb = colorbar;
clb.Ticks =linspace(0.0,1, 10);
clb.TickLabels = T_(1:2:end).*1000;
clb.Label.String = '$\tau$ (ms)';
clb.Direction = 'reverse';
clb.Label.Interpreter = 'latex';
clb.TickLabelInterpreter = 'latex';
clb.FontSize = 12;
clb.Label.FontSize = 14;
colormap(cMap2(1:length(T_), :));

%% print figure
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Simulations_frequency','-dpdf','-r0')
%% plot filgure 4
T__ = linspace(0,1,3000).';
Freqs = (1:20)'*1000;
CovFreq = exp( - magicFactor * (Freqs .* sqrt(T__).'*volatility).^2 );

f = figure(6);clf; hold on; grid on; box on

plot(1000*T__,CovFreq.^2);
colororder(cMap2)
% plot(T,Cor(100:100:400,:));
xlabel('Time (ms)', 'Interpreter','latex')
ylabel('$R^2_\vartheta (\tau, f)$', 'Interpreter','latex')
set(gca, 'FontSize', 12)
clb = colorbar;
clb.Ticks = [1000; Freqs(2:2:end)]./20000 -0.025;
clb.TickLabels = [1000; Freqs(2:2:end)]./1000;
clb.Label.String = 'Frequency (kHz)';
clb.Direction = 'reverse';
clb.Label.Interpreter = 'latex';
clb.TickLabelInterpreter = 'latex';
clb.FontSize = 12;
clb.Label.FontSize = 14;
colormap(cMap2);

%% print figure 
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Simulations_time','-dpdf','-r0')

%% Simulate averaging process

avg_s = mean(s,3);

avgEnergy = abs(avg_s).^2;
avgEnergy = movsum(avgEnergy,freqWin,1);
avgEnergy = movsum(avgEnergy,timeWin,2);

avgEnergyLoss_dB = pow2db(avgEnergy ./ Energy);

avgEnergyLossModel_db = pow2db((CovModel+1)/2); 

figure;
surf(T,F,avgEnergyLoss_dB,'EdgeColor','none')
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
view([0 90])
colorbar

figure; hold on; grid on;
plot(T,avgEnergyLossModel_db(100:100:400,:));
set(gca,'ColorOrderIndex',1);
plot(T,avgEnergyLoss_dB(100:100:400,:));
xlabel('Time (seconds)')
ylabel('Energy loss (dB)')

%% synthesize fractional delays
function rir = synthRIR(echo_times, echo_amps, fs, len)

lenSamples = ceil(len*fs);
numberOfEchoes = size(echo_times,1);
RIR = zeros(lenSamples,1);

impulse = zeros(lenSamples,1);
impulse(2) = 1;

Impulse = fft(impulse);
A = 1i*angle(Impulse);

for it = 1:numberOfEchoes
    RIR = RIR + echo_amps(it) .* exp(A.* (echo_times(it)*fs-1));
end

rir = real(ifft(RIR));

end