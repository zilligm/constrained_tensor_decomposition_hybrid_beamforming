%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HYBRID BEAMFORMING SIMULATIONS
%   Author: Guilherme Martignago Zilli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Simulation Parameters
    NR      = 500;                   % Number of realizations         
    SNR_dB  = -20 : 5 : 20;         % SNR values in dB
    SNR     = 10.^(SNR_dB./10);     % SNR values in decimal
        
 % Scenario Parameters
        K       = 1;                % # of users (K = 1 for SU, K>1 for MU) 
        Ns      = 4;                % # or data stream (per user)
        M       = 512;               % # of subcarriers (for OFDM only)  - NOT IMPLEMENTED!

    % Transmitter side
        Nt      = 64;               % # of antennas at transmitter
        Nt_rf   = K*Ns;             % # of RF chains at transmitter
    
    % Receiver side
        Nr      = 64;               % # of antennas at receiver
        Nr_rf   = Ns;               % # of RF chains at receiver
    
    % Channel Parameters
        Ncl     = 5;                % # of clusters
        Nray    = 10;               % # of rays in each cluster
  
    % Total Transmitted Power
        Pt = K*Ns;
      
% Checking conditions:
    if ((K*Ns > Nt_rf) || (Nt_rf > Nt))
        error('Check conditions at BS')
    end
    if ((Ns > Nr_rf) || (Nr_rf > Nr))
        error('Check conditions at MS')
    end
    
% Codebook for Chiang2018
%     Nf = 10; 
    CDBK = kron(dftMatrix(sqrt(Nt),'linear'),dftMatrix(sqrt(Nt),'linear'));

    
% Parallel Variable Initialization
    % Spectral Efficiency
        Ropt         = zeros(length(SNR),1);
        Rsoh         = zeros(length(SNR),1);
        Rtsai        = zeros(length(SNR),1);
        Rpe          = zeros(length(SNR),1);
        Rch          = zeros(length(SNR),1);
        Rjzhang      = zeros(length(SNR),1);
        Rmy          = zeros(length(SNR),1);

    % Execution time
        Tsoh         = zeros(1,1);
        Ttsai        = zeros(1,1);
        Tpe          = zeros(1,1);
        Tch          = zeros(1,1);
        Tmy          = zeros(1,1);
    
% Algorithm Simulations
parfor nreal = 1:NR
    nreal
    
    % Spectral Efficiency
        ropt         = zeros(length(SNR),1);
        rsoh         = zeros(length(SNR),1);
        rtsai        = zeros(length(SNR),1);
        rpe          = zeros(length(SNR),1);
        rch          = zeros(length(SNR),1);
        rjzhang      = zeros(length(SNR),1);
        rmy          = zeros(length(SNR),1);
        
    % Execution time
        tsoh         = zeros(length(SNR),1);
        ttsai        = zeros(1,1);
        tpe          = zeros(1,1);
        tch          = zeros(length(SNR),1);
        tmy          = zeros(1,1);

        F_pe = zeros(Nt,Ns,M);
        W_pe = zeros(Nr,Ns,M);
        F_my = zeros(Nt,Ns,M);
        W_my = zeros(Nr,Ns,M);
        

    % Channel Matrix
    [H,At,Ar,Fopt,Wopt] = channel_realization(Nt,Nr,K,M,Ns,Ncl,Nray,'square',Pt);
        
    % Precoder & Combiner Design    

        % Sohrabi2017
            % Algorithms runs at every SNR
            
        % THTsai 2019
            tstart = tic;
            [ F_tsai, W_tsai ] = THTsai2019( H, Nt_rf, Nr_rf, Ns, Pt);
            ttsai = toc(tstart);
            
        % Phase Extraction [Yu2016]   
            tstart = tic;
            [FRF_pe, FBB_pe] = PE_AltMin(Fopt, Nt_rf);
                for m = 1:M
                    FBB_pe(:,:,m) = sqrt(Pt) * FBB_pe(:,:,m) / norm(FRF_pe * FBB_pe(:,:,m),'fro');
                    F_pe(:,:,m) = FRF_pe*FBB_pe(:,:,m);
                end
            [WRF_pe, WBB_pe] = PE_AltMin(Wopt, Nr_rf);
                for m = 1:M
                    W_pe(:,:,m) = WRF_pe*WBB_pe(:,:,m);
                end
            tpe = toc(tstart);
            
        % MyBeamforming
            tstart = tic;
                [ F_my, W_my, DATA ] = myBeamformingFastTucker( H, Ns, Pt);
            tmy = toc(tstart);
            
        % JZhang2016
            [ F_jzhang, W_jzhang] = JZhang2016( H, Ns, Pt );                  

                
        for snr_ind = 1:length(SNR)
            sigma2 = 0.000000001;
            rho = sigma2*SNR(snr_ind)/Pt;
            
            %Sohrabi2017
                tstart = tic;
                    [F_Soh,W_Soh] = Sohrabi2017(H,Ns, sigma2,Pt);
                tsoh = toc(tstart);  
                
            % Chiang2018
                tstart = tic;
                    [ F_ch, W_ch] = Chiang2018( H, Ns, CDBK, Ns, sigma2, Pt );                
                tch(snr_ind) = toc(tstart);
                
            % Spectral Efficiency
                ropt(snr_ind)         = SUspectralEfficiency(H,OFDMWaterFilling(H,Fopt,Wopt,rho,'total',Pt),Wopt,Ns,rho,sigma2);
                rsoh(snr_ind)         = SUspectralEfficiency(H,OFDMWaterFilling(H,F_Soh,W_Soh,rho,'total',Pt),W_Soh,Ns,rho,sigma2);
                rtsai(snr_ind)        = SUspectralEfficiency(H,OFDMWaterFilling(H,F_tsai,W_tsai,rho,'total',Pt),W_tsai,Ns,rho,sigma2);
                rpe(snr_ind)          = SUspectralEfficiency(H,OFDMWaterFilling(H,F_pe,W_pe,rho,'total',Pt),W_pe,Ns,rho,sigma2);
                rch(snr_ind)          = SUspectralEfficiency(H,OFDMWaterFilling(H,F_ch,W_ch,rho,'total',Pt),W_ch,Ns,rho,sigma2);
                rjzhang(snr_ind)      = SUspectralEfficiency(H,OFDMWaterFilling(H,F_jzhang,W_jzhang,rho,'total',Pt),W_jzhang,Ns,rho,sigma2);
                rmy(snr_ind)          = SUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);   
        end
    
    % Spectral Efficiency
        Ropt         = Ropt + ropt;
        Rsoh         = Rsoh + rsoh;
        Rtsai        = Rtsai + rtsai;
        Rpe          = Rpe + rpe;
        Rch          = Rch + rch;
        Rjzhang      = Rjzhang + rjzhang;
        Rmy          = Rmy + rmy;
        
    % Execution time
        Tsoh         = Tsoh + tsoh/NR;
        Ttsai        = Ttsai + ttsai/NR;
        Tpe          = Tpe + tpe/NR;
        Tch          = Tch + max(tch)/NR;
        Tmy          = Tmy + tmy/NR;
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

% filename = 'Figures\OFDM_SE_vs_SNR_data';
% save(filename)


% Spectral Efficiency
h = figure
hold on; grid on
plot(SNR_dB,Ropt,'-ok','LineWidth',1.5);
plot(SNR_dB,Rsoh,'-ob','LineWidth',1.5);
plot(SNR_dB,Rtsai,'-oy','LineWidth',1.5);
plot(SNR_dB,Rpe,'-oc','LineWidth',1.5);
plot(SNR_dB,Rch,'-om','LineWidth',1.5);
plot(SNR_dB,Rjzhang,'-og','LineWidth',1.5);
plot(SNR_dB,Rmy,'-or','LineWidth',1.5);
legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)')
xlabel('SNR (dB)')

% filenameFigFile = 'Figures\OFDM_SE_vs_SNR_data';
% saveas(gcf,filenameFigFile,'pdf');