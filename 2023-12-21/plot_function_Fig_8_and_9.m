%% plot
figure; hold on; grid on;
plot(bandCenters,divergence')
xlabel('Frequency Band (Hz)');
ylabel('Divergence (seconds)');

% plots only first two RIR correlation
figure; hold on; grid on;
offset = 0*(1:size(pred_cor,2));
plot( time_cor, squeeze(meas_cor(:,1,:)) + offset)
set(gca,'ColorOrderIndex',1);
plot( pred_toastd ./ (sqrt(2)*divergence(1,:)), pred_cor + offset)
xlim([0 max(time_cor)])


%% colors
numPlots = 5;

colorMod = linspace(1,0,numPlots);
col1 = [0, 0.4470, 0.7410];
cMap = [col1(1) * colorMod; col1(2)*colorMod; col1(3)*colorMod];
cred = [1 0 0];
cVec1 = linspace(0,1, numPlots);
cMap2 = [cVec1; col1(2)*colorMod; col1(3)*colorMod];

col2 = [113, 62, 90]./255;

%%

for it = 1:5
eNoise(:, it) = (energy(0.7*length(energy):0.8*length(energy),it));

end
eN = median(eNoise, 1);

for it = 1:5
    snr(:, it) = (energy(:, it)-eN(it))/eN(it);
end
snrdB = db(snr);
[~, locs] = min(abs(snrdB(1: 20000, :)));

exp_corr  = snr./(1+snr);

%% plot HV
figure(1); clf; hold on
ind = 19;
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
for i = 1:4;
    plot( time_cor*1000, squeeze(meas_cor(:,i,ind)) , 'color', cMap2(:, i), 'LineWidth',1)
plot( 1000*pred_toastd ./ (sqrt(2)*divergence(i,ind)), pred_cor(:, ind) ,'--', 'color', cMap2(:, i), 'HandleVisibility','off', 'LineWidth',1)
end
plot(1000*[mean(locs)./fs, mean(locs)./fs], [-.4 1], 'k--')
% ylim([0.9 1.005])
xlim([0 0.4].*1000);
xlabel('Time (ms)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation', 'Interpreter','latex', 'FontSize',12)
% set(gca, 'ytick', .9:.02:1, 'xtick', 0:.04:.2, 'FontSize',12)
% lgd  = legend('$1$', '$2$', '$3$', '$4$', 'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'southwest', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 2);
lgd.Title.String = [ {'Time between'},{' measurements'}];%$j-i$';
box on
grid on

%% get the correlation coefficient between the calculated correlation curves and the fitted curves
clear new_cor

cor_PCC = zeros(2,2,4,bandIt);
for i = 1:4
    for it = 1: bandIt
        firstI = find(mask(:, it) == 1, 1, 'first');
        firstI = firstI+500;
        corL = find(mask(firstI+1:end, it) == 0, 1, 'first');

        pred_time(:, i, it) = pred_toastd ./ (sqrt(2)*divergence(i,it));
        [~, loc] = min(abs(pred_time(1:end, i, it) - corL/fs));
        new_cor = interp1(pred_time(1:loc, i, it), pred_cor(1:loc,  it), time_cor(firstI:corL));
        new_cor = new_cor(~isnan(new_cor));
        nc_len = length(new_cor);
        cor_PCC (:, :, i, it) = corrcoef([meas_cor(firstI:firstI+nc_len-1,i, it), new_cor]);
    end
end
%% plot the correlation coefficient for the fits
f = figure(7); clf; hold on

for i =1:4 
    plot(bandCenters, squeeze(cor_PCC(1,2,i, :)), '-.o', 'color', cMap2(:, i), 'LineWidth',1, 'MarkerFaceColor',cMap2(:, i))
end 
yline(0.7, 'k--', 'HandleVisibility','off')
% ylim([0.75 1])
 xlim(1000*[0.5 19.5]);
xlabel('Center frequency (kHz)', 'Interpreter','latex', 'FontSize',12)
ylabel('Correlation coefficient', 'Interpreter','latex', 'FontSize',12)
set(gca,  'xtick', 1000*[1:2:20], 'xticklabel', 1:2:20,  'FontSize',12)
% set(gca, 'ytick',[0.75:0.1:1, 1])
set(gca, 'ytick',[-0.5 0 0.5 0.7 1])
box on
f.Position(end) = 250;

lgd  = legend('5 s', '10 s', '15 s', '20 s', 'location', 'southeast', 'interpreter', 'latex', 'fontsize', 12, 'numcolumns', 4);
lgd.Title.String = [ {'Time between measurements'}];%$j-i$';
%% print the correlation coefficient figure
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[f.Position(3), f.Position(4)])
print(f,'corr_coef','-dpdf','-r0')