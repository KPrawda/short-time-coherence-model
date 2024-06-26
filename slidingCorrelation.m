function [r_corr, e_sig, r_snr, e_ref, e_mean] = slidingCorrelation(irRef, ir, winLen)
%% calculate the RIR coherence 
% % Input:
% % - irRef = reference RIR
% % - ir = RIRs to calculate coherence
% % - winLen = length of the analysis window in samples
% % Output:
% % - r_corr = calculated coherence
% % - e_sig = RIR energy
% % - r_snr = SNR-based expected coherence
% % - e_ref = energy of the reference RIR
% % - e_mean = energy of the averaged RIRs
num_rirs = size(ir,2);

% estimate noise level
noise = irRef(round(0.9*end):end,:);
noiseLevel = sqrt(mean(noise.^2,1));

% compute correlation
mov = @(x) movsum(x,winLen)./winLen;

% win = hann(winLen);
win = rectwin(winLen);
% win = blackman(winLen);
win = win / sum(win);
mov = @(x) conv(x,win,'same');

for it = 1:num_rirs
    r_cov(:,it) = mov(irRef.*ir(:,it)); % covariance
    e_sig(:,it) = mov(ir(:,it).^2);     % signal energy

    sig = [irRef, ir(:,it)];
    avg_s = mean(sig,2);                % average the signal
    avgEnergy = abs(avg_s).^2;          % energy of the averaged signal
    avgEnergy = mov(avgEnergy);         
    e_mean(:, it) =  avgEnergy;
end
e_ref =  mov(irRef.^2);                 % energy of the reference RIR
r_corr = r_cov./sqrt(e_sig.*e_ref);     % coherence between a RIR and a reference RIR

% estimate rir energy
e_rir = max(e_sig - noiseLevel.^2 * 1,0);

% SNR-based correlation
r_snr = e_rir ./ e_sig;


end