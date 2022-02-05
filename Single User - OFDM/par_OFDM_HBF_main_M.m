%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HYBRID BEAMFORMING SIMULATIONS
%   Author: Guilherme Martignago Zilli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all
clc

% Simulation Parameters
    NR      = 500;                 % Number of realizations  
    SNR_dB  = 5;                    % SNR values in dB
    SNR     = 10.^(SNR_dB./10);     % SNR values in decimal

    % Scenario Parameters
        MM      = 2.^[2:9]; 
        K       = 1;                            % # of users (K = 1 for SU, K>1 for MU) \
        Ns      = 4;            % # of data stream (per user)
    % Transmitter side
        Nt      = 8^2;              % # of antennas at transmitter
        Nt_rf   = Ns;                % # of RF chains at transmitter
    % Receiver side
        Nr      = 8^2;               % # of antennas at receiver
        Nr_rf   = Ns;                % # of RF chains at receiver
    % Channel Parameters
        Ncl     = 5;                % # of clusters
        Nray    = 10;               % # of rays in each cluster

        Pt = K*Ns;
        sigma2 = 0.000000001;
        rho = sigma2*SNR/Pt; 
        
% Codebook for Chiang2018
%     Nf = 10; 
    CDBK = kron(dftMatrix(sqrt(Nt),'linear'),dftMatrix(sqrt(Nt),'linear'));
        
% Parallel Variable Initialization
    % Spectral Efficiency
    Ropt         = zeros(length(MM),1);
    Rsoh         = zeros(length(MM),1);
    Rtsai        = zeros(length(MM),1);
    Rpe          = zeros(length(MM),1);
    Rch          = zeros(length(MM),1);
    Rjzhang      = zeros(length(MM),1);
    Rmy          = zeros(length(MM),1);
        
% Algorithm Simulations

parfor nr = 1:NR
    nr
    % Spectral Efficiency
        ropt         = zeros(length(MM),1);
        rsoh         = zeros(length(MM),1);
        rtsai        = zeros(length(MM),1);
        rpe          = zeros(length(MM),1);
        rch          = zeros(length(MM),1);
        rjzhang      = zeros(length(MM),1);
        rmy          = zeros(length(MM),1);
        
    for s = 1:length(MM)
        M = MM(s);

        % Checking conditions:
        if ((K*Ns > Nt_rf) || (Nt_rf > Nt))
            error('Check conditions at BS')
        end
        if ((Ns > Nr_rf) || (Nr_rf > Nr))
            error('Check conditions at MS')
        end
        
        F_pe = zeros(Nt,Ns,M);
        W_pe = zeros(Nr,Ns,M);
        F_my = zeros(Nt,Ns,M);
        W_my = zeros(Nr,Ns,M);    

        % Channel Matrix
            [H,At,Ar,Fopt,Wopt] = channel_realization(Nt,Nr,K,M,Ns,Ncl,Nray,'square',Pt);

        % Sohrabi2017
            [F_Soh,W_Soh] = Sohrabi2017(H,Ns,sigma2,Pt);

        % THTsai 2019
            [ F_tsai, W_tsai ] = THTsai2019( H, Nt_rf, Nr_rf, Ns, Pt);

        % Phase Extraction [Yu2016]
            [FRF_pe, FBB_pe] = PE_AltMin(Fopt, Nt_rf);
                for m = 1:M
                    FBB_pe(:,:,m) = sqrt(Pt) * FBB_pe(:,:,m) / norm(FRF_pe * FBB_pe(:,:,m),'fro');
                    F_pe(:,:,m) = FRF_pe*FBB_pe(:,:,m);
                end
            [WRF_pe, WBB_pe] = PE_AltMin(Wopt, Nr_rf);
                for m = 1:M
                    W_pe(:,:,m) = WRF_pe*WBB_pe(:,:,m);
                end

        % Chiang2018
            [ F_ch, W_ch] = Chiang2018( H, Ns, CDBK, Ns, sigma2, Pt ); 

        % JZhang2016
            [ F_jzhang, W_jzhang] = JZhang2016( H, Ns, Pt );

        % MyBeamforming
            [ F_my, W_my, DATA ] = myBeamformingFastTucker( H, Ns, Pt );      

        % Spectral Efficiency
            ropt(s)         = SUspectralEfficiency(H,OFDMWaterFilling(H,Fopt,Wopt,rho,'total',Pt),Wopt,Ns,rho,sigma2);
            rsoh(s)         = SUspectralEfficiency(H,OFDMWaterFilling(H,F_Soh,W_Soh,rho,'total',Pt),W_Soh,Ns,rho,sigma2);
            rtsai(s)        = SUspectralEfficiency(H,OFDMWaterFilling(H,F_tsai,W_tsai,rho,'total',Pt),W_tsai,Ns,rho,sigma2);
            rpe(s)          = SUspectralEfficiency(H,OFDMWaterFilling(H,F_pe,W_pe,rho,'total',Pt),W_pe,Ns,rho,sigma2);
            rch(s)          = SUspectralEfficiency(H,OFDMWaterFilling(H,F_ch,W_ch,rho,'total',Pt),W_ch,Ns,rho,sigma2);
            rjzhang(s)      = SUspectralEfficiency(H,OFDMWaterFilling(H,F_jzhang,W_jzhang,rho,'total',Pt),W_jzhang,Ns,rho,sigma2);
            rmy(s)          = SUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);     
    end
    
    % Spectral Efficiency
        Ropt         = Ropt + ropt;
        Rsoh         = Rsoh + rsoh;
        Rtsai        = Rtsai + rtsai;
        Rpe          = Rpe + rpe;
        Rch          = Rch + rch;
        Rjzhang      = Rjzhang + rjzhang;
        Rmy          = Rmy + rmy;

end

% Spectral Efficiency
Ropt         = Ropt/NR;
Rsoh         = Rsoh/NR;
Rtsai        = Rtsai/NR;
Rpe          = Rpe/NR;
Rch          = Rch/NR;
Rjzhang      = Rjzhang/NR;
Rmy          = Rmy/NR;

%% Results
legendCell{1} = 'Fully digital (SVD)';
legendCell{2} = '[Sohrabi2017]';
legendCell{3} = 'SS-SVD [THTsai2019]';
legendCell{4} = 'PE-AltMin [Yu2016]';
legendCell{5} = 'ICSI-HBF [Chiang2018]';
legendCell{6} = 'PE-HOSVD [JZhang2016]';
legendCell{7} = 'Proposed HBF';

% filename = strcat('Figures\OFDM_SE_vs_M_data_SNR_',num2str(SNR_dB));
% save(filename)

h = figure
plot(MM,Ropt,'-ok','LineWidth',1.5); hold on; grid on
plot(MM,Rsoh,'-ob','LineWidth',1.5);
plot(MM,Rtsai,'-oy','LineWidth',1.5);
plot(MM,Rpe,'-oc','LineWidth',1.5);
plot(MM,Rch,'-om','LineWidth',1.5);
plot(MM,Rjzhang,'-og','LineWidth',1.5);
plot(MM,Rmy,'-or','LineWidth',1.5);
legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)')
xlabel('Number of subcarrier - $M$')
xlim([min(MM) max(MM)])

% filenameFigFile = strcat('Figures\OFDM_SE_vs_M_fig_SNR_',num2str(SNR_dB));
% saveas(gcf,filenameFigFile,'pdf');