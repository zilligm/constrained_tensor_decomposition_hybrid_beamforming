clear all
clc

% Simulation Parameters
    NR      = 1000;                 % Number of realizations  
%     SNR_dB  = -20 : 5 : 10;          % SNR values in dB
%     SNR     = 10.^(SNR_dB./10);     % SNR values in decimal

% Scenario Parameters
        K       = 8;                % # of users (K = 1 for SU, K>1 for MU) 
        Ns      = 8;                % # or data stream (per user)
        % Nf      = 1;                % # of subcarriers (for OFDM only)  - NOT IMPLEMENTED!

    % Transmitter side
        Nt      = 144;              % # of antennas at transmitter
        Nt_rf   = K*Ns;                % # of RF chains at transmitter
    
    % Receiver side
        Nr      = 16;               % # of antennas at receiver
        Nr_rf   = Ns;                % # of RF chains at receiver
    
    % Channel Parameters
        Ncl     = 5;                % # of clusters
        Nray    = 10;               % # of rays in each cluster
        
    % Array Type
        arrayType = 'square';
        
        if ( strcmp(arrayType,'square') && ~( (~( mod(Nt,sqrt(Nt)) )) || (~( mod(Nr,sqrt(Nr)) )) ) )
            error('The number of antennas must be N = sqrt(N)*sqrt(N) with sqrt(N) integer when using USPA')
        end

% Checking conditions:
    if ((K*Ns > Nt_rf) || (Nt_rf > Nt))
        error('Check conditions at BS')
    end
    if ((Ns > Nr_rf) || (Nr_rf > Nr))
        error('Check conditions at MS')
    end
    
    
    NITE    = zeros(Ns,NR);
            
% Algorithm Simulations
for nreal = 1:NR
    nreal
%     [H,At,Ar,Fopt,Wopt,D] = channel_realization(Nt,Nr,K,Ns,Ncl,Nray,arrayType);
    H = randn(Nr,Nt,K) + 1i*randn(Nr,Nt,K);
    [nite] = myMUBeamformingConv( H, Ns );    
    NITE(:,nreal) = nite;             
end
MinMin = min(min(NITE))
MaxMax = max(max(NITE))
MeanMin = mean(min(NITE))
MeanMax = mean(max(NITE))
MeanMean = mean(mean(NITE,1))

%% Results
[fmin,xmin] = ecdf(min(NITE)); 
[fmax,xmax] = ecdf(max(NITE));
[favg,xavg] = ecdf(mean(NITE)); 


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


 