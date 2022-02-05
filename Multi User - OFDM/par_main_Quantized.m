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
    SNR_dB  = -20 : 5 : 20;          % SNR values in dB
    SNR     = 10.^(SNR_dB./10);     % SNR values in decimal

% Scenario Parameters
        K       = 4;                % # of users (K = 1 for SU, K>1 for MU) 
        Ns      = 4;                % # or data stream (per user)
        M       = 32;
        Pt      = K*Ns;

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
    
    sinf   = zeros(length(SNR),1);
    s6     = zeros(length(SNR),1);
    s5     = zeros(length(SNR),1);
    s4     = zeros(length(SNR),1);
    s3     = zeros(length(SNR),1);
    s2     = zeros(length(SNR),1);
    s1     = zeros(length(SNR),1); 
 
    % Channel Matrix
    [H,At,Ar,Fopt,Wopt,D] = channel_realization(Nt,Nr,K,M,Ncl,Nray,arrayType);

    % Precoder & Combiner Design
    for snr_ind = 1:length(SNR)
        sigma2 = 0.000000001;
        rho = sigma2*SNR(snr_ind)/Pt;

        % Fully Digital SVD-based
%             [F_svd,W_svd] = uncontrainedSVD_MUBF( H, Ns );  

        % My MU Beamforming
            [F_inf,W_inf] = myMUOFDMBeamformingRCD(H,Ns, sigma2, Pt, rho); 
            
        % 6-bits quantization
            [F_6b,W_6b] = myMUOFDMBeamformingRCD_Quant(H,Ns, sigma2, Pt, rho,6); 
        
        % 5-bits quantization
            [F_5b,W_5b] = myMUOFDMBeamformingRCD_Quant(H,Ns, sigma2, Pt, rho,5); 
        
        % 4-bits quantization
            [F_4b,W_4b] = myMUOFDMBeamformingRCD_Quant(H,Ns, sigma2, Pt, rho,4); 
        
        % 3-bits quantization
            [F_3b,W_3b] = myMUOFDMBeamformingRCD_Quant(H,Ns, sigma2, Pt, rho,3); 
            
        % 2-bits quantization
            [F_2b,W_2b] = myMUOFDMBeamformingRCD_Quant(H,Ns, sigma2, Pt, rho,2);  
        
        % 1-bits quantization
            [F_1b,W_1b] = myMUOFDMBeamformingRCD_Quant(H,Ns, sigma2, Pt, rho,1); 

        % Sum Sum-Rate with WaterFilling
            sinf(snr_ind)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_inf,W_inf,rho,'total',Pt),W_inf,Ns,rho,sigma2);
            s6(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_6b,W_6b,rho,'total',Pt),W_6b,Ns,rho,sigma2);
            s5(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_5b,W_5b,rho,'total',Pt),W_5b,Ns,rho,sigma2);
            s4(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_4b,W_4b,rho,'total',Pt),W_4b,Ns,rho,sigma2);
            s3(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_3b,W_3b,rho,'total',Pt),W_3b,Ns,rho,sigma2);
            s2(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_2b,W_2b,rho,'total',Pt),W_2b,Ns,rho,sigma2);
            s1(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_1b,W_1b,rho,'total',Pt),W_1b,Ns,rho,sigma2);
    end
    
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

% filename = 'Figures\MU_SR_vs_SNR_Quant_data';
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

% filenameFigFile = 'Figures\MU_SE_vs_SNR_Quant_fig';
% saveas(gcf,filenameFigFile,'pdf');