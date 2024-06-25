%% plot room with colored speed of sound voxels - Figure 1 from the manuscript 
%% and the resulting RIRs - Figure 2 from the manuscript
% complementary code for the publication 
% "Short-time Coherence Between Repeated Room Impulse Response Measurements"
% by K. Prawda, S. J. Schlecht, and V. Välimäki
% submitted to the Journal of the Acoustical Society of America
% on 22.03.2024
%% housekeeping
clear all; close all; clc
%% colors
numPlots = 10;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;
%% set x and y axes
x=0:1:5;
y = 0:1:3;
[X,Y] = meshgrid(x,y);
%% plot Fig. 1a
rng(1) % for reproducibility
T = 273.15 + 20 + 3*(rand(size(X))-0.5); % vary temperature by +- 1.5 deg C
C = 331.45 + sqrt(1 + T./273); % calculate speed of sound

f = figure(1); clf
subplot(1,2,1)

hold on
surf(X, Y, C, FaceAlpha=.9, FaceColor="flat", EdgeColor="none", EdgeAlpha=0)% speed of sound voxels in the real room
surf(fliplr(X-5), Y, C, FaceAlpha=.7, FaceColor="flat", EdgeColor="none", EdgeAlpha=0) % speed of sound voxels in the image room 1
surf((X), flipud(Y+3), C, FaceAlpha=.7, FaceColor="flat", EdgeColor="none", EdgeAlpha=0) % speed of sound voxels in the image room 2

% plot walls
plot3([0 5],[0 0],[1 1], 'k', 'linewidth', 2)
plot3([0 5],[3 3],[1 1], 'k', 'linewidth', 2)
plot3([0 0],[0 3],[1 1], 'k', 'linewidth', 2)
plot3([5 5],[0 3],[1 1], 'k', 'linewidth', 2)

%Receiver
plot3(4,1, 400, 'wo' , 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','w') % plot receiver 
text(4.2 ,1, 1.1, 'R', 'color', 'w', 'Interpreter','latex', 'FontSize',14)% letter R by the receiver

%Real source
scatter3(1,2, 400,150, 'k','filled', 'o', 'markerfacealpha', 1)% plot real sound source
text(0.9,2.3, 400, 'S', 'color', 'w', 'Interpreter','latex', 'FontSize', 14) % letter S by the sound source

%Image source 1
scatter3(-1,2, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1, 'MarkerEdgeColor',[.25 .25 .25])
text(-1.0,2.3, 400, 'IS$_1$', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)

%Image source 2
scatter3(1,4, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1, 'MarkerEdgeColor',[.25 .25 .25])
text(1,4.3, 400, 'IS$_2$', 'color', 'w', 'Interpreter','latex', 'FontSize',14)

%sound paths 
plot3([1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
plot3([-1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
plot3([1 4], [4 1], 400*[1 1], 'k--', 'linewidth',1)

% time label
text(-5 ,5.75, 1.1, '$c(\textit{\textbf{x}})$', 'color', 'k', 'Interpreter','latex', 'FontSize',14)

axis tight
x0=100;
y0=100;
width=0.7*1000;
height=0.7*600;
set(gcf,'position',[x0,y0,width,height])
set(gca, 'ytick', [], 'xtick', [])
grid off
box off
view(2)
axis off

% color option and colorbar
colormap(cMap2.')
clb =  colorbar('location', 'southoutside');
clim([332.884 332.896])
clb.TickLabelInterpreter = 'latex';
clb.Ticks = [332.885 332.890 332.895];
clb.Label.String = [{'Speed of sound (m/s)'},{'(a)'}];
clb.Label.Interpreter = 'Latex';
clb.FontSize = 14;

% get rid of the axes
set(gca,'ycolor','w', 'xcolor','w')

%% plot Fig. 1b
rng(10) % for reproducibility
C_ =C +  0.005*(2*rand(size(X))-1); % vary the speed of sound by a random value

subplot(1,2,2); 
hold on
% plot room
surf(X, Y, C_, FaceAlpha=.9,  EdgeColor="none", EdgeAlpha=0) % speed of sound voxels in the real room
surf(fliplr(X-5), Y, C_, FaceAlpha=.7, EdgeColor="none", EdgeAlpha=0) % speed of sound voxels in the image room 1
surf((X), flipud(Y+3), C_, FaceAlpha=.7,  EdgeColor="none", EdgeAlpha=0) % speed of sound voxels in the image room 2

% plot walls
plot3([0 5],[0 0],[1 1], 'k', 'linewidth', 2)
plot3([0 5],[3 3],[1 1], 'k', 'linewidth', 2)
plot3([0 0],[0 3],[1 1], 'k', 'linewidth', 2)
plot3([5 5],[0 3],[1 1], 'k', 'linewidth', 2)

% Receiver
plot3(4,1, 400, 'wo' , 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','w')
text(4.2 ,1, 1.1, 'R', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)

% Real source
scatter3(1,2, 400,150, 'k','filled', 'o', 'markerfacealpha', 1)
text(0.9,2.3, 400, 'S', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)

% Image source 1
scatter3(-1,2, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1, 'MarkerEdgeColor',[.25 .25 .25])
text(-1.0,2.3, 400, 'IS$_1$', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)

% Image source 2
scatter3(1,4, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1, 'MarkerEdgeColor',[.25 .25 .25])
text(1,4.3, 400, 'IS$_2$', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)

% Sound paths
plot3([1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
plot3([-1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
plot3([1 4], [4 1], 400*[1 1], 'k--', 'linewidth',1)

% Time label
text(-5 ,5.75, 1.1, '$c^{\prime}(\textit{\textbf{x}})$', 'color', 'k', 'Interpreter','latex', 'FontSize',14)


% axis and figure 
axis tight
x0=100;
y0=100;
width=1000*0.7;
height=600*0.7;
set(gcf,'position',[x0,y0,width,height])
set(gca, 'ytick', [], 'xtick', [])
grid off
box off
view(2)
axis off

% colormap and colorbar
colormap(cMap2.')
clb =  colorbar('location', 'southoutside');
clim([332.884 332.896])
clb.TickLabelInterpreter = 'latex';
clb.Ticks = [332.885 332.890 332.895];
clb.Label.String = [{'Speed of sound (m/s)'},{'(b)'}];
clb.Label.FontSize = 14;
clb.Label.Interpreter = 'Latex';
clb.FontSize =  14;



set(gca,'ycolor','w', 'xcolor','w')
f.Position = [100 100 1500 420];
%% print the ISM figure (Fig. 1)
h=gcf;
set(h,'PaperPositionMode','auto');
set(h,'PaperOrientation','landscape');
set(h,'Position',f.Position);
exportgraphics(h,'ISM_1_2.jpg','Resolution',600)
%% load the pre-computed ISM
load("IR_ISM.mat")

time = 0:1/fs:length(h.ism)/fs - 1/fs;
%% plot Fig 2 from the paper 
f = figure(2); clf;  hold on; box on; 
plot(1000*time, h.ism(:, 1), 'color',cMap2(:,1), 'LineWidth',2) % plot time in ms
plot(1000*time, h.ism(:, 2), 'color',cMap2(:,end),  'LineWidth',1)
ylim([-0.025 0.075])
xlim([0 400])
f.Position(end) = 220;
rectangle('Position', [200 -0.003 5 0.01],  'LineWidth',1)

% plot the dashed lines from rectangle in IR to the inlay axes
plot([200 165], [0.007 0.06], 'k--',  'LineWidth',1)
plot([200 165], [-0.003 0.021], 'k--',  'LineWidth',1)
plot([205 370], [0.007 0.06], 'k--',  'LineWidth',1)
plot([205 370], [-0.003 0.021], 'k--',  'LineWidth',1)

% 
xlabel('Time (ms)', 'Interpreter','latex', 'fontsize', 12)
ylabel('Signal value', 'Interpreter','latex', 'fontsize', 12);
set(gca, 'fontsize', 12)
legend('RIR at $c(\textit{\textbf{x}})$', 'RIR at $c^{\prime}(\textit{\textbf{x}})$',  'Interpreter','latex', 'fontsize', 12, 'location', 'southeast', 'numcolumns',2)
% inlay axis with zoom on the IR
ax = axes('position',[.45 .55 .4 .3]);
plot(1000*time, h.ism(:, 1), 'color',cMap2(:,1), 'LineWidth',2); hold on
plot(1000*time, h.ism(:, 2), 'color',cMap2(:,end),  'LineWidth',1.5)
xlim([200.2 201])
ylim([-0.002 0.006])
% box on
ax.YTick = [];
ax.XTick = [];
%% print the IR figure
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
% print(f,'IR_difference','-dpdf','-r0')