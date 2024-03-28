%% complementary code for publication 
% Scḧafer, M., Prawda, K., Rabenstein, R., Schlecht, S.J., "Distribution
% of Modal Damping in Absorptive Shoebox Rooms", in Proc. WASPAA 2023, New Paltz, NY, USA, Oct 22--25, 2023

% Compare proposed method to ISM and Lehmann methods
% and plot Figure 5
% uses:
% * shoebox2modes() 
% * histwv() 
% * modalSynthesis() 
% * zeroPhaseBandpass()
% * compute_echograms_mics()    (in './shoebox-roomsim')
% * render_mic_rirs()           (in './shoebox-roomsim')
% * rir_generator()             (in './RIR-Generator/')
% * fast_ISM_RoomResp()         (in './fastISM')
% * edc()

% Sebastian J. Schlecht, Thursday, 13 April 2023
% shoebox to damping density comparison to ISM
%% housekeeping 
clear; clc; close all;
addpath('C:\Users\prawdak1\Dropbox (Aalto)\Projects\Loudspeaker deconvolution\ISM scripts/shoebox-roomsim')
addpath('C:\Users\prawdak1\Dropbox (Aalto)\Projects\Loudspeaker deconvolution\ISM scripts/RIR-Generator/')

%% set simulation parameters
fs = 48000;          % Sample frequency (samples/s)

% room dimensions
Lx = 5;   % in [m]
Ly = 3;   % in [m]
Lz = 3;   % in [m]
L = [Lx, Ly, Lz];

c = [332.888; 332.892];    % speed of sound in [m/s]

rec = [4 1 1];            % Receiver position [x y z] (m)
src = [1 2 1];            % Source position [x y z] (m)
limitsTime = 1.2;               % lenght of impulse response (seconds)
time = (1:limitsTime*fs).'/fs;  % seconds

bandpassEdges = [50 20000];     % Hz

numberOfCases = 2;%3;

%% simulation
for it = 1:numberOfCases
    % reflection coefficients
    a = 0.2*ones(1,6);
    r = sqrt(1-a);
    

    % image sources
    switch 'Habets'
        case 'Politis'
            [abs_echograms, rec_echograms, echograms] = compute_echograms_mics(L, src, rec, a, limitsTime);
            h_temp = render_mic_rirs(abs_echograms, [], fs);
            h.ism(:,it) = h_temp(1:limitsTime*fs,:);
        case 'Habets'
            mtype = 'omnidirectional';    % Type of microphone
            order = -1;                 % -1 equals maximum reflection order!
            dim = 3;                    % Room dimension
            orientation = 0;     % Microphone orientation (rad)
            hp_filter = 0;              % Disable high-pass filter
            h.ism(:,it) = rir_generator(c(it), fs, rec, src, L, r, limitsTime*fs, mtype, order, dim, orientation, hp_filter).';
    end
    

       
end

h.ism = zeroPhaseBandpass(h.ism,bandpassEdges,fs);


%% colors for plotting
c2 =[239, 45, 86]./255;
c1 = [46, 134, 171]./255;
c3 = [255, 203, 119]./255;
%% plot Figure 5
close all;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
cuton = 0;%0.2*fs;
f = figure(1); hold on; grid on;

plot(time,h.ism./max(abs(h.ism)),'-', 'LineWidth',1);


set(gca, 'FontSize',12)

xlabel('Time (seconds)', 'Interpreter','latex', 'FontSize',12)
ylabel('Signal value', 'Interpreter','latex', 'FontSize',12);

box on
% legend(lines, {'ISM'}, 'Interpreter','latex', 'FontSize',12)
f.Position(end) = 270;

%%
save('IR_ISM.mat', 'h', 'fs')