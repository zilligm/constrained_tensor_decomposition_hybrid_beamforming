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
        M       = 512;               % # of subcarriers (for OFDM only) 

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
    
    Sinf   = zeros(length(SNR),1);
    S6     = zeros(length(SNR),1);
    S5     = zeros(length(SNR),1);
    S4     = zeros(length(SNR),1);
    S3     = zeros(length(SNR),1);
    S2     = zeros(length(SNR),1);
    S1     = zeros(length(SNR),1);

    
% Algorithm Simulations
parfor nreal = 1:NR
    nreal
    
    % Spectral Efficiency
        sinf   = zeros(length(SNR),1);
        s6     = zeros(length(SNR),1);
        s5     = zeros(length(SNR),1);
        s4     = zeros(length(SNR),1);
        s3     = zeros(length(SNR),1);
        s2     = zeros(length(SNR),1);
        s1     = zeros(length(SNR),1);
       
    % Channel Matrix
    [H,~,~,~,~] = channel_realization(Nt,Nr,K,M,Ns,Ncl,Nray,'square',Pt);
                
        % My MU Beamforming
            [F_inf,W_inf,~] = myBeamformingFastTucker( H, Ns, Pt);
            
        % 6-bits quantization
            [F_6b,W_6b,~] = myBeamformingFastTucker_Quant( H, Ns, 6, Pt);
        
        % 5-bits quantization
            [F_5b,W_5b,~] = myBeamformingFastTucker_Quant( H, Ns, 5, Pt);
        
        % 4-bits quantization
            [F_4b,W_4b,~] = myBeamformingFastTucker_Quant( H, Ns, 4, Pt);
        
        % 3-bits quantization
            [F_3b,W_3b,~] = myBeamformingFastTucker_Quant( H, Ns, 3, Pt);
            
        % 2-bits quantization
            [F_2b,W_2b,~] = myBeamformingFastTucker_Quant( H, Ns, 2, Pt);
        
        % 1-bits quantization
            [F_1b,W_1b,~] = myBeamformingFastTucker_Quant( H, Ns, 1, Pt);
      
        for snr_ind = 1:length(SNR) 
            sigma2 = 0.000000001;
            rho = sigma2*SNR(snr_ind)/Pt;
            
            % Sum Sum-Rate with WaterFilling
            sinf(snr_ind)   = SUspectralEfficiency(H,OFDMWaterFilling(H,F_inf,W_inf,rho,'total',Pt),W_inf,Ns,rho,sigma2);
            s6(snr_ind)     = SUspectralEfficiency(H,OFDMWaterFilling(H,F_6b,W_6b,rho,'total',Pt),W_6b,Ns,rho,sigma2);
            s5(snr_ind)     = SUspectralEfficiency(H,OFDMWaterFilling(H,F_5b,W_5b,rho,'total',Pt),W_5b,Ns,rho,sigma2);
            s4(snr_ind)     = SUspectralEfficiency(H,OFDMWaterFilling(H,F_4b,W_4b,rho,'total',Pt),W_4b,Ns,rho,sigma2);
            s3(snr_ind)     = SUspectralEfficiency(H,OFDMWaterFilling(H,F_3b,W_3b,rho,'total',Pt),W_3b,Ns,rho,sigma2);
            s2(snr_ind)     = SUspectralEfficiency(H,OFDMWaterFilling(H,F_2b,W_2b,rho,'total',Pt),W_2b,Ns,rho,sigma2);
            s1(snr_ind)     = SUspectralEfficiency(H,OFDMWaterFilling(H,F_1b,W_1b,rho,'total',Pt),W_1b,Ns,rho,sigma2);
        end
    
    % Spectrak Efficiency 
    Sinf   = Sinf + sinf;
    S6     = S6 + s6;
    S5     = S5 + s5;
    S4     = S4 + s4;
    S3     = S3 + s3;
    S2     = S2 + s2;
    S1     = S1 + s1;
end 

Sinf   = Sinf/NR;
S6     = S6/NR;
S5     = S5/NR;
S4     = S4/NR;
S3     = S3/NR;
S2     = S2/NR;
S1     = S1/NR;

%% Results
legendCell{1} = 'Nb = inf';
legendCell{2} = 'Nb = 6';
legendCell{3} = 'Nb = 5';
legendCell{4} = 'Nb = 4';
legendCell{5} = 'Nb = 3';
legendCell{6} = 'Nb = 2';
legendCell{7} = 'Nb = 1';

% filename = 'Figures\SU_SR_vs_SNR_Quant_data';
% save(filename);


figure
hold on; grid on
plot(SNR_dB,Sinf,'r-o','LineWidth',1.5);
plot(SNR_dB,S6,'r-x','LineWidth',1.5);
plot(SNR_dB,S5,'c-o','LineWidth',1.5);
plot(SNR_dB,S4,'g-o','LineWidth',1.5);
plot(SNR_dB,S3,'b-o','LineWidth',1.5);
plot(SNR_dB,S2,'k-o','LineWidth',1.5);
plot(SNR_dB,S1,'m-o','LineWidth',1.5);

legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)')
xlabel('SNR (dB)')

% filenameFigFile = 'Figures\SU_SE_vs_SNR_Quant_fig';
% saveas(gcf,filenameFigFile,'pdf');