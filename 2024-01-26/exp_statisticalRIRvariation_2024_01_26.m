% Sebastian J. Schlecht, Tuesday, 16. January 2024
% test time-frequency correlation model based on normal TDoA distribution
clear; clc; close all;

fs = 48000*1;
lenSeconds = 1;
lenSamples = fs*lenSeconds;

% generate echo times
num_echoes = 2*10^4;
mean_echo_times = rand(num_echoes,1) * lenSeconds;

% generate variation
divergence = 10^-6*3; % TOASTD increase per second.
num_rirs = 2;

toaVariation = sqrt(mean_echo_times) .* divergence .* randn(num_echoes,num_rirs); % Gaussian random walk - is the divergance correct here or should it be squared?
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

%% the model of time-frequency correlation is based on the Gaussian TDoA

magicFactor = 20; % might be influenced by the STFT parameters
CovModel = exp( - magicFactor * (F .* sqrt(T).'*divergence).^2 );

figure;
surf(T,F,CovModel,'EdgeColor','none')
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
view([0 90])
colorbar

figure; hold on; grid on;
plot(T,CovModel(100:100:400,:));
set(gca,'ColorOrderIndex',1);
plot(T,Cor(100:100:400,:));
xlabel('Time (seconds)')
ylabel('Correlation')

figure; hold on; grid on;
plot(F,CovModel(:,300:100:500));
set(gca,'ColorOrderIndex',1);
plot(F,Cor(:,300:100:500));
xlabel('Frequency (Hz)')
ylabel('Correlation')

%% Simulate averaging process

avg_s = mean(s,3);

avgEnergy = abs(avg_s).^2;
avgEnergy = movsum(avgEnergy,freqWin,1);
avgEnergy = movsum(avgEnergy,timeWin,2);

avgEnergyLoss_dB = pow2db(avgEnergy ./ Energy);

avgEnergyLossModel_db = pow2db((CovModel+1)/2); % formula for the function of energy loss - this is for compemsation model later

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