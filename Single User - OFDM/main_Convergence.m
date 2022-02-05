%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HYBRID BEAMFORMING SIMULATIONS
%   Author: Guilherme Martignago Zilli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Simulation Parameters
    NR      = 1000;                 % Number of realizations  

% Scenario Parameters
        K       = 1;                % # of users (K = 1 for SU, K>1 for MU) 
        Ns      = 8;                % # or data stream (per user)
        M       = 64;                % # of subcarriers (for OFDM only)  - NOT IMPLEMENTED!

    % Transmitter side
        Nt      = 64;               % # of antennas at transmitter
        Nt_rf   = K*Ns;             % # of RF chains at transmitter
    
    % Receiver side
        Nr      = 64;               % # of antennas at receiver
        Nr_rf   = Ns;               % # of RF chains at receiver
    
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
    
% Algorithm Simulations
TESTE_MY = zeros(Ns,NR);

for nr = 1:NR
    nr

    % Channel Matrix
    [H,At,Ar,Fopt,Wopt] = channel_realization(Nt,Nr,K,M,Ns,Ncl,Nray,'linear');
        % MyBeamforming
            [ FRF_my, WRF_my, DATA ] = myBeamforming( H, Ns );
            
            TESTE_MY(:,nr) = DATA;
                
end    
%% Results
[fmin,xmin] = ecdf(min(TESTE_MY)); 
[fmax,xmax] = ecdf(max(TESTE_MY));
[favg,xavg] = ecdf(mean(TESTE_MY)); 

min(min(TESTE_MY))
mean(mean(TESTE_MY))    
max(max(TESTE_MY))
%% Results
figure
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
h = legend('min(\eta)','average(\eta)','max(\eta)','Interpreter','Latex','Location','best')
h.Position = [0.6   0.35    0.25    0.25];
h.FontName = 'Times New Roman';
h.FontSize = 8.5;

