%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HYBRID BEAMFORMING SIMULATIONS
%   Author: Guilherme Martignago Zilli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Simulation Parameters
    NR      = 500;                 % Number of realizations  
%     SNR_dB  = -20 : 5 : 10;          % SNR values in dB
%     SNR     = 10.^(SNR_dB./10);     % SNR values in decimal

% Scenario Parameters
        K       = 1;                % # of users (K = 1 for SU, K>1 for MU) 
        Ns      = 4;                % # or data stream (per user)
        M       = 512;                % # of subcarriers (for OFDM only)  - NOT IMPLEMENTED!

    % Transmitter side
        Nt      = 64;              % # of antennas at transmitter
        Nt_rf   = K*Ns;                % # of RF chains at transmitter
    
    % Receiver side
        Nr      = 64;               % # of antennas at receiver
        Nr_rf   = Ns;                % # of RF chains at receiver
    
    % Channel Parameters
        Ncl     = 5;                % # of clusters
        Nray    = 10;               % # of rays in each cluster
    

% Checking conditions:
    if ((K*Ns > Nt_rf) || (Nt_rf > Nt))
        error('Check conditions at BS')
    end
    if ((Ns > Nr_rf) || (Nr_rf > Nr))
        error('Check conditions at MS')
    end
    
    DATA = zeros(Ns,NR);
    
% Algorithm Simulations
for nreal = 1:NR
    nreal

    % Channel Matrix
    [H,At,Ar,Fopt,Wopt] = channel_realization(Nt,Nr,K,M,Ns,Ncl,Nray,'square',K*Ns);

    [ F_my, W_my, data ] = myBeamformingFastTucker( H, Ns, K*Ns );

    DATA(:,nreal) = data;

end   
[fmin,xmin] = ecdf(min(DATA)); 
[fmax,xmax] = ecdf(max(DATA));
[favg,xavg] = ecdf(mean(DATA)); 

minDATA = min(min(DATA))
meanDATA = mean(mean(DATA))    
maxDATA = max(max(DATA))

clear DATA
clear F_my
clear W_my
clear At Ar Fopt Wopt H
%% Results
[minDATA meanDATA maxDATA]

h = figure
fig_width = 9;
fig_heigth = 5;
    fig_pos = [4, 2, fig_width, fig_heigth]; % Figura com tamanho 4.5x9cm
    fig_config = {'Color', [1,1,1], 'Units', 'Centimeters','PaperType', 'A4',...
                    'PaperOrientation','portrait','Position', fig_pos,...
                    'PaperPositionMode', 'auto',};
    fig=gcf;
    set(fig,fig_config{:});
    legend(get(fig,'CurrentAxes'),'off');

ax = gca;
ax.ColorOrderIndex = 1;
plot(xmin,fmin,'-.','LineWidth',1.5,'Color','b'); hold on
plot(xavg,favg,':','LineWidth',1.5,'Color','k');
plot(xmax,fmax,'--','LineWidth',1.5,'Color','r');
ylabel('Empirical CDF')
xlabel('Number of Iterations')
fig.CurrentAxes.FontName = 'Times New Roman';
% fig.CurrentAxes.XLim = [0 40];
fig.CurrentAxes.Position = [0.119    0.175    0.82    0.76];
% h = legend('min(N_{ite})','mean(N_{ite})','max(N_{ite})','Interpreter','Latex','Location','best');
leg = legend('min(\eta)','average(\eta)','max(\eta)','Interpreter','Latex','Location','best')
leg.Position = [0.6   0.35    0.25    0.25];
leg.FontName = 'Times New Roman';
leg.FontSize = 8.5;


% filename = 'Figures\OFDM_convergence_Nmax=100';
% save(filename)
% filenameFigFile = filename;
% saveas(h,filenameFigFile,'fig')
 
% filename = 'Figures\OFDM_convergence';
% save(filename)