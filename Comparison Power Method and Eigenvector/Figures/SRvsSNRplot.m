clc
clear all
close all

data = open('PM_vs_EG_data.mat')
close all

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesFontSize',9); 
set(groot, 'defaultLineLineWidth',1.3); 
set(groot, 'defaultLineMarkerSize',7);

figure;
setFig(6.5,9.75);
fig = gcf;
ax = gca;
fig.CurrentAxes.Position = [0.13    0.15    0.82    0.82];
ax = gca;

ax.ColorOrderIndex = 1;
plot(data.SNR_dB,data.SRmyPM,'-o','LineWidth',1.5); hold on; grid on
plot(data.SNR_dB,data.SRmyEV,'-+','LineWidth',1.5);
ylabel('Sum-Rate(bits/s/Hz)','Interpreter','latex')
xlabel('SNR (dB)','Interpreter','latex')
% xlim([data.SNR_dB(2) data.SNR_dB(end-1)])
xlim([data.SNR_dB(1) data.SNR_dB(end)])
ylim([0 140])
yticks([0:20:140])

h = legend('Power Method','Eigenvector','Location','best','Interpreter','latex')
% h.Position = [0.13   0.8468    0.3432    0.1238];
h.Position = [0.130    0.8466    0.3187    0.1235];


% create a new pair of axes inside current figure
    axes('position',[0.7775    0.3399    0.1500    0.2000])
%     axes('position',[ 0.2376    0.4253    0.1500    0.2000])
    ax = gca;
    ax.ColorOrderIndex = 1;
    box on % put box around new pair of axes
        %ax.ColorOrderIndex = 1;
        plot(data.SNR_dB,data.SRmyPM,'-o','LineWidth',1.5); hold on; grid on
        plot(data.SNR_dB,data.SRmyEV,'-+','LineWidth',1.5);
    fig.CurrentAxes.YLim = [56.6 57.4];
    fig.CurrentAxes.XLim = [-0.06 0.06];
    fig.CurrentAxes.XAxis.TickValues = 0;
    fig.CurrentAxes.YAxis.TickValues = fig.CurrentAxes.YLim;
% 

save = 1;
if save==1
    print(fig,'-depsc',strcat('Figures\PowerMethVsEigenvector_SR','.eps'));
    print(fig,'-dpdf',strcat('Figures\PowerMethVsEigenvector_SR','.pdf'));
end    
% close





function setFig(fig_heigth,fig_width)
   
    fig_pos = [4, 2, fig_width, fig_heigth]; % Figura com tamanho 4.5x9cm
    fig_config = {'Color', [1,1,1], 'Units', 'Centimeters','PaperType', 'A4',...
                    'PaperOrientation','portrait','Position', fig_pos,...
                    'PaperPositionMode', 'auto',};
    fig=gcf;
    set(fig,fig_config{:});
    legend(get(fig,'CurrentAxes'),'off');

end


