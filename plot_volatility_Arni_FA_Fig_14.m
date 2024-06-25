%% plotting Arni 2 volatility - Figure 14 from the manuscript
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024
%% housekeeping
clear all
close all
clc
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
%% load pre-calculated volatility data
load('Arni_volatility_FA.mat');
%%
bandCenters = 1:19; %kHz
numMeas = max(size(volatility)); 
%% create colormap
numPlots =10;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cMap3 = [linspace(col1(1),1,numPlots); linspace(col1(2),1,numPlots); linspace(col1(3),1,numPlots)];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;

%% median filtering
med_div = medfilt1(volatility, 50,[], 1, "omitnan", "truncate");
%% plot the figure 
f = figure(1); clf; hold on

% subplot(2,1,1); hold on
for it = 10:18%bandCenters
    scatter([1:numMeas]*5/60, squeeze(volatility( :, it)),10,cMap2(:,it-10+1)','o', 'filled', 'markerfacealpha', .07,'MarkerEdgeColor',cMap2(:,it-10+1)','MarkerEdgeAlpha', 0.12,'handlevisibility', 'off')
    plot([1:numMeas]*5/60, med_div(:, it),'color', cMap2(:,it-10+1)','linewidth', 2.5)
    
end

ylim([1.5*10^-6, 0.8*10^-5])
xlim([-50 numMeas*5+50]./60)
box on
set(gca, 'fontsize', 12)
xlabel('Time between measurements (min)', 'Interpreter','latex', 'fontsize', 12)
ylabel('Volatility $\vartheta\,(\textrm{s}/\sqrt{\textrm{s}})$', 'Interpreter','latex', 'fontsize', 12);

colormap(cMap2(:, 1:end-1)')
caxis([10 18])
clb = colorbar;
clb.Label.String = 'Frequency (kHz)';
clb.Label.Interpreter = 'latex';
clb.Ticks = linspace(10.5, 17.5,9);
clb.TickLabels = 10:18;
clb.TickLabelInterpreter = 'latex';
% clb.Limits = [1, 19];
clb.FontSize = 12;
f.Position(end) = 350;


%% print figure
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Arni_FA_divergence','-dpdf','-r0')