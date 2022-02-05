clc
clear all
close all

data = open('PM_vs_EG_data.mat')
MINrE = min( min(vec(data.NITEpm(:,:,end))) , min(vec(data.NITEev(:,:,end))) );
MAXrE = max( max(vec(data.NITEpm(:,:,end))) , max(vec(data.NITEev(:,:,end))) );

rangE = [MINrE-0.5:1:MAXrE+0.5];
close all

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesFontSize',9); 
set(groot, 'defaultLineLineWidth',1.3); 
set(groot, 'defaultLineMarkerSize',7);

figure;
setFig(6.5,10.5);
fig = gcf;
ax = gca;
fig.CurrentAxes.Position = [0.13    0.15    0.75    0.82];
ax = gca;


% yyaxis left
% ax.ColorOrderIndex = 1;
% histogram(vec(data.NITE),data.rangE,'Normalization','probability')
% ylabel('Empirical PMF','Interpreter','latex');
% xlabel('Number of iterations','Interpreter','latex');
% xlim([0 60])
% 
% yyaxis right
% ax.ColorOrderIndex = 2;
% [favgall,xavgall] = ecdf(vec(data.NITE));
% plot(xavgall,favgall,':','LineWidth',1.5);
% ylabel('Empirical CDF','Interpreter','latex');
% % xlabel('Number of Iterations','Interpreter','latex');
% % xlim([0 100])


[favgallNB,xavgallNB] = ecdf(vec(data.NITEpm(:,:,:)));
[favgallOFDM,xavgallOFDM] = ecdf(vec(data.NITEev(:,:,:)));
plot(xavgallNB,favgallNB,':'); hold on
plot(xavgallOFDM,favgallOFDM,':'); 
ylabel('Empirical CDF','Interpreter','latex');
xlabel('Number of Iterations','Interpreter','latex');
xlim([0 60])


h = legend('Power Method','Eigenvector','Location','best','Interpreter','latex')
% h.Position = [0.13   0.8468    0.3432    0.1238];
h.Position = [0.5594    0.1493    0.3187    0.1235];

grid on




save = 1;
if save==1
    print(fig,'-depsc','-r600',strcat('Figures\PowerMethVsEigenvectorCDF','.eps'));
    print(fig,'-dpdf','-r600',strcat('Figures\PowerMethVsEigenvectorCDF','.pdf'));
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


