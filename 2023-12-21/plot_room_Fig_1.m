%% plot room 
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
%%
x=0:1:5
y = 0:1:3;

[X,Y] = meshgrid(x,y);


%%
% C = C./max(abs(C));
rng(1)
T = 273.15 + 20 + 3*(rand(size(X))-0.5);
C = 331.45 + sqrt(1 + T./273);%2*rand(size(X))-1;%C +rand(size(A))-0.5;
f = figure(1); clf
subplot(1,2,1)
% set(gcf, 'color', 'none');   
% subplot(2,2,4); hold on
% s = pcolor(X, Y, C);
% s.FaceAlpha=.9, s.FaceColor="flat", s.EdgeColor="none", s.EdgeAlpha=0
hold on
surf(X, Y, C, FaceAlpha=.9, FaceColor="flat", EdgeColor="none", EdgeAlpha=0)
surf(fliplr(X-5), Y, C, FaceAlpha=.7, FaceColor="flat", EdgeColor="none", EdgeAlpha=0)


surf((X), flipud(Y+3), C, FaceAlpha=.7, FaceColor="flat", EdgeColor="none", EdgeAlpha=0)

plot3([0 5],[0 0],[1 1], 'k', 'linewidth', 2)
plot3([0 5],[3 3],[1 1], 'k', 'linewidth', 2)
plot3([0 0],[0 3],[1 1], 'k', 'linewidth', 2)
plot3([5 5],[0 3],[1 1], 'k', 'linewidth', 2)

plot3(4,1, 400, 'wo' , 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','w')


scatter3(1,2, 400,150, 'k','filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
text(0.9,2.3, 400, 'S', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)
text(-1.0,2.3, 400, 'IS$_1$', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)
scatter3(-1,2, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])

scatter3(1,4, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% % plot3(1,4, 1, 'o' ,'color', [.4 .4 .4 ], 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(1,4, 1,'ko', 'linewidth',1, 'markersize',12)
% % scatter3(1.2,1.95, 1,150, 'w','filled', 'o', 'markerfacealpha', .15, 'markeredgecolor', 'none')% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(1.2,1.95, 1,'ko', 'linestyle', ':', 'markersize',12, 'linewidth',1)
plot3([1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
% scatter3(-1,2, 1,150,  [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(-1,2, 1,'ko', 'linewidth',1, 'markersize',12)
% % scatter3(-1.1,2.02, 1,150, 'w','filled', 'o', 'markerfacealpha', .1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(-1.1,2.02, 1,'ko', 'linestyle', '--', 'markersize',12, 'linewidth',1)
plot3([-1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
% % scatter3(0.8,4.2, 1,150, 'w','filled', 'o', 'markerfacealpha', .1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ]
% plot3(0.8,4.2, 400,'ko', 'linestyle', '--', 'markersize',12, 'linewidth',1)
text(1,4.3, 400, 'IS$_2$', 'color', 'w', 'Interpreter','latex', 'FontSize',14)
plot3([1 4], [4 1], 400*[1 1], 'k--', 'linewidth',1)
% scatter3(1,4, 1,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(4,1, 1.1, 'wo' , 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','w')
text(4.2 ,1, 1.1, 'R', 'color', 'w', 'Interpreter','latex', 'FontSize',14)
text(-5 ,5.75, 1.1, 'Time $t_1$', 'color', 'k', 'Interpreter','latex', 'FontSize',14)

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

colormap(cMap2.')
clb =  colorbar('location', 'southoutside');
caxis([332.888 332.892])
clb.TickLabelInterpreter = 'latex'
clb.Ticks = [332.888 332.890 332.892]
% clb.TickLabels = {'Low', 'Average', 'High'}
clb.Label.String = [{'Speed of sound'},{'(a)'}]
clb.Label.Interpreter = 'Latex'
% set(clb, 'Direction', 'reverse')
clb.FontSize = 14

% subplot(2,2,3)


set(gca,'ycolor','w', 'xcolor','w')

%%
% C = C./max(abs(C));
rng(10)
C_ =C +  0.002*(2*rand(size(X))-1);%C +rand(size(A))-0.5;
% f = figure(2); clf
subplot(1,2,2); %clf
% set(gcf, 'color', 'none');   
% subplot(2,2,4); hold on
hold on
surf(X, Y, C_, FaceAlpha=.9,  EdgeColor="none", EdgeAlpha=0)
surf(fliplr(X-5), Y, C_, FaceAlpha=.7, EdgeColor="none", EdgeAlpha=0)
surf((X), flipud(Y+3), C_, FaceAlpha=.7,  EdgeColor="none", EdgeAlpha=0)

plot3([0 5],[0 0],[1 1], 'k', 'linewidth', 2)
plot3([0 5],[3 3],[1 1], 'k', 'linewidth', 2)
plot3([0 0],[0 3],[1 1], 'k', 'linewidth', 2)
plot3([5 5],[0 3],[1 1], 'k', 'linewidth', 2)

plot3(4,1, 400, 'wo' , 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','w')


scatter3(1,2, 400,150, 'k','filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
text(0.9,2.3, 400, 'S', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)
scatter3(-1,2, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
text(-1.0,2.3, 400, 'IS$_1$', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)
scatter3(1,4, 400,150, [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
text(1,4.3, 400, 'IS$_2$', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)
plot3(1,4, 1,'ko', 'linewidth',1, 'markersize',12)
% plot3(1,4, 1, 'o' ,'color', [.4 .4 .4 ], 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])

% scatter3(1.3,1.92, 1,150, 'w','filled', 'o', 'markerfacealpha', .25, 'markeredgecolor', 'none')% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(1.3,1.92, 1,'ko', 'linestyle', '--', 'markersize',12, 'linewidth',1)
plot3([1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)

% scatter3(-1,2, 1,150,  [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% scatter3(-1.2,2.05, 1,150, 'w','filled', 'o', 'markerfacealpha', .25)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(-1.2,2.05, 1,'ko', 'linestyle', '--', 'markersize',12, 'linewidth',1)
% plot3(-1,2, 1,'ko', 'markersize',12, 'linewidth',1)
plot3([-1 4], [2 1], 400*[1 1], 'k--', 'linewidth',1)
% scatter3(0.95,4.03, 1,150, 'w','filled', 'o', 'markerfacealpha', .25)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ]
% plot3(0.95,4.03, 1,'ko', 'linestyle', '--', 'markersize',12, 'linewidth',1)
plot3([1 4], [4 1], 400*[1 1], 'k--', 'linewidth',1)
% scatter3(1,4, 1,150,  [.5 .5 .5],'filled', 'o', 'markerfacealpha', 1)% 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor', [.2 .2 .2 ])
% plot3(4,1, 1.1, 'wo' , 'LineWidth',2, 'MarkerSize',10, 'MarkerFaceColor','w')
text(4.2 ,1, 1.1, 'R', 'color', 'w', 'Interpreter','latex', 'FontSize', 14)
text(-5 ,5.75, 1.1, 'Time $t_2$', 'color', 'k', 'Interpreter','latex', 'FontSize',14)


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
colormap(cMap2.')
clb =  colorbar('location', 'southoutside');
caxis([332.888 332.892])
clb.TickLabelInterpreter = 'latex'
% clb.Ticks = [-1 0 1]
clb.Ticks = [332.888 332.890 332.892]
% clb.TickLabels = {'Low', 'Average', 'High '}
clb.Label.String = [{'Speed of sound'},{'(b)'}]
clb.Label.FontSize = 14
clb.Label.Interpreter = 'Latex'
% set(clb, 'Direction', 'reverse')
clb.FontSize =  14
% caxis([-1 1])

% subplot(2,2,3)


set(gca,'ycolor','w', 'xcolor','w')
f.Position = [100 100 1500 420]
%%
h=gcf;
set(h,'PaperPositionMode','auto');
set(h,'PaperOrientation','landscape');
set(h,'Position',f.Position);
exportgraphics(h,'ISM_1_2.jpg','Resolution',600)
%% 
load("IR_ISM.mat")

%% upsample the IRs for better visibility
upFactor = 1;
h_up = upsample(h.ism, upFactor);
time = 0:1/(upFactor*fs):(length(h.ism))/fs - 1/(upFactor*fs);


%%
f = figure(2); clf;  hold on; box on; 
plot(1000*time, h_up(:, 1), 'color',cMap2(:,1), 'LineWidth',2)
plot(1000*time, h_up(:, 2), 'color',cMap2(:,end),  'LineWidth',1)
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
legend('RIR at $t_1$', 'RIR at $t_2$',  'Interpreter','latex', 'fontsize', 12, 'location', 'southeast', 'numcolumns',2)
% inlay axis
ax = axes('position',[.45 .55 .4 .3]);
plot(1000*time, h_up(:, 1), 'color',cMap2(:,1), 'LineWidth',2); hold on
plot(1000*time, h_up(:, 2), 'color',cMap2(:,end),  'LineWidth',1)
xlim([200 200.5])
ylim([-0.002 0.006])
% box on
ax.YTick = [];
ax.XTick = [];
%% print the IR figure
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'IR_difference','-dpdf','-r0')