%% plotting Arni divergence
clear all; close all; clc

%%
load('Arni_volatility.mat')
%% colors
numPlots = 5;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;

%% plot
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

%%
v_s = volatility;
center_freqs = 1:19;
%%
n1 = [];
n2 = [];
n3 = [];
n4 = [];

f = figure(6);clf; hold on;
for i = 1:5
    rirInd = 1:5;
    rirInd(i) = [];
    for j = abs(rirInd-i)
            scatter((center_freqs+(j-2)/5)',squeeze(v_s(:, i, j, :)),20,cMap2(:, 5-j)','o', 'filled', 'markerfacealpha', .07,'handlevisibility', 'off') %,'color', cMap2(:, i), 'LineWidth',1.5) ;
            if j ==1
                n1 = [n1; squeeze(v_s(:, i, j, :))];
            elseif j == 2
                n2 = [n2; squeeze(v_s(:, i, j, :))];
            elseif j == 3
                n3 = [n3; squeeze(v_s(:, i, j, :))];
            elseif j == 4
                n4 = [n4; squeeze(v_s(:, i, j, :))];
            
            end         
    end
end
%%
scatter(center_freqs( :) -1/5, median(n1, 1),32,cMap2(:,  4).','o', 'filled', 'markerfacealpha', 1, 'markeredgecolor', [1 1 1], 'linewidth', 0.5)
scatter(center_freqs( :)-0/5 , median(n2, 1),32,cMap2(:,  3).','o', 'filled', 'markerfacealpha', 1, 'markeredgecolor', [1 1 1], 'linewidth', 0.5) 
scatter(center_freqs( :)+1/5, median(n3, 1),32,cMap2(:,  2).','o', 'filled', 'markerfacealpha', 1, 'markeredgecolor', [1 1 1], 'linewidth', 0.5) 
scatter(center_freqs( :) +2/5, median(n4, 1),32,cMap2(:,  1).','o', 'filled', 'markerfacealpha', 1, 'markeredgecolor', [1 1 1], 'linewidth', 0.5) 
ylim([1*10^-7, 1*10^-5])
set(gca, 'xtick', [1:2:20], 'fontsize', 12)
labs = [0.2:0.2:1.2].*10^(-5);
xlim([0.2 19.5])
xlabel('Center frequency (kHz)', 'Interpreter','latex', 'fontsize', 12)
ylabel('Volatility $(s/\sqrt{s})$', 'Interpreter','latex', 'fontsize', 12);
box on

lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'northwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 4);
lgd.Title.String = [ {'Time between measurements'}];


%% print figure 
set(f,'Units','Inches');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'Arni_divergence','-dpdf','-r0')

