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
        M       = 32;                % # of subcarriers (for OFDM only)  - NOT IMPLEMENTED!
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
    
    SmyPM      = zeros(length(SNR),1);
    SmyEV      = zeros(length(SNR),1);
    
    SRmyPM      = zeros(length(SNR),1);
    SRmyEV      = zeros(length(SNR),1);
    
    NITEpm = zeros(K*Ns,NR);
    NITEev = zeros(K*Ns,NR);
    
% Algorithm Simulations
for nreal = 1:NR
    nreal

    smyPM      = zeros(length(SNR),1);    
    smyEV      = zeros(length(SNR),1);    

    srmyPM      = zeros(length(SNR),1);
    srmyEV      = zeros(length(SNR),1);    
    
    
    % Channel Matrix
    [H,~,~,~,~,~] = channel_realization(Nt,Nr,K,M,Ncl,Nray,arrayType); 

    for snr_ind = 1:length(SNR)        
        rho= 1;
        sigma2 = rho*Pt/SNR(snr_ind);
        
        [F_my,W_my] = myMUOFDMBeamformingRCD(H,Ns, sigma2,Pt);
            smyPM(snr_ind)       = MUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);
            srmyPM(snr_ind)       = MUsumRate(H,OFDMWaterFilling(H,F_my,W_my,SNR(snr_ind),'total'),W_my,Ns,rho,sigma2);
        
        [F_my,W_my] = myMUOFDMBeamformingRCD_eig(H,Ns, sigma2,Pt); 
            smyEV(snr_ind)       = MUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);
            srmyEV(snr_ind)       = MUsumRate(H,OFDMWaterFilling(H,F_my,W_my,SNR(snr_ind),'total'),W_my,Ns,rho,sigma2);
    end
    NITEpm(:,nreal,snr_ind) = DATAofdm;
    NITEev(:,nreal,snr_ind) = DATAnarrow;
       
    Sfdbd    = Sfdbd + sfdbd;
    Sgonz    = Sgonz + sgonz;
    Sdzhang  = Sdzhang + sdzhang;
    Smy      = Smy + smy;
    
    SRfdbd    = SRfdbd + srfdbd;
    SRgonz    = SRgonz + srgonz;
    SRdzhang  = SRdzhang + srdzhang;
    SRmy      = SRmy + srmy;
    
end   
Sfdbd    = Sfdbd/NR;
Sgonz    = Sgonz/NR;
Sdzhang  = Sdzhang/NR;
Smy      = Smy/NR;

SRfdbd    = SRfdbd/NR;
SRgonz    = SRgonz/NR;
SRdzhang  = SRdzhang/NR;
SRmy      = SRmy/NR;

%% Results
legendCell{1} = 'FD-BD';
legendCell{end+1} = 'Gonzalez2018';
legendCell{end+1} = 'DZhang2019';
legendCell{end+1} = 'Prop. Alg.';

filename = 'Figures\MU_OFDM_SR_vs_SNR_data';
save(filename);


figure
hold on; grid on
plot(SNR_dB,Sfdbd,'k-o','LineWidth',1.5);
plot(SNR_dB,Sgonz,'g-x','LineWidth',1.5);
plot(SNR_dB,Sdzhang,'b-o','LineWidth',1.5);
plot(SNR_dB,Smy,'r-o','LineWidth',1.5);

legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)','interpreter','latex')
xlabel('SNR (dB)','interpreter','latex')

filenameFigFile = 'Figures\MU_OFDM_SE_vs_SNR_fig';
saveas(gcf,filenameFigFile,'pdf');

figure
hold on; grid on
plot(SNR_dB,SRfdbd,'k-o','LineWidth',1.5);
plot(SNR_dB,SRgonz,'g-x','LineWidth',1.5);
plot(SNR_dB,SRdzhang,'b-o','LineWidth',1.5);
plot(SNR_dB,SRmy,'r-o','LineWidth',1.5);

legend(legendCell,'Location','northwest')
ylabel('Average Sum-Rate (bits/s/Hz)','interpreter','latex')
xlabel('SNR (dB)','interpreter','latex')

filenameFigFile = 'Figures\MU_OFDM_SR_vs_SNR_fig';
saveas(gcf,filenameFigFile,'pdf');