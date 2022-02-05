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
    
    Sfdbd    = zeros(length(SNR),1);
    Sfdrbd    = zeros(length(SNR),1);
    Sgonz    = zeros(length(SNR),1);
    Sdzhang  = zeros(length(SNR),1);
    Smy      = zeros(length(SNR),1);
    
    
% Algorithm Simulations
parfor nreal = 1:NR
    nreal
    
    sfdbd    = zeros(length(SNR),1);
    sfdrbd    = zeros(length(SNR),1);
    sgonz    = zeros(length(SNR),1);
    sdzhang  = zeros(length(SNR),1);
    smy      = zeros(length(SNR),1);    
    
    % Channel Matrix
    [H,~,~,~,~,~] = channel_realization(Nt,Nr,K,M,Ncl,Nray,arrayType);
    
    % Fully Digital Block Diagonalization
        [F_fdbd,W_fdbd] = FDBD(H,Ns,Pt);
        
    for snr_ind = 1:length(SNR)
        sigma2 = 0.000000001;
        rho = sigma2*SNR(snr_ind)/Pt;

        [F_fdrbd,W_fdrbd] = FDRBD(H,Ns,Pt,sigma2,rho);
        
%         [ F_Gonz, W_Gonz] = GonzalesComa2018D1(H, Ns, sigma2 ,Pt, rho);
%            
        [ F_DZhang, W_DZhang] = DZhang2019_TensorUnfolding( H, Ns, sigma2 ,Pt, rho);

        [F_my,W_my] = myMUOFDMBeamformingRCD(H,Ns, sigma2,Pt, rho); 
             
        % Spectral Efficiency with WaterFilling
            sfdbd(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_fdbd,W_fdbd,rho,'total',Pt),W_fdbd,Ns,rho,sigma2);
            sfdrbd(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_fdrbd,W_fdrbd,rho,'total',Pt),W_fdrbd,Ns,rho,sigma2);
%             sgonz(snr_ind)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_Gonz,W_Gonz,rho,'total',Pt),W_Gonz,Ns,rho,sigma2);
            sdzhang(snr_ind)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_DZhang,W_DZhang,rho,'total',Pt),W_DZhang,Ns,rho,sigma2);
            smy(snr_ind)       = MUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);
    end

       
    Sfdbd    = Sfdbd + sfdbd;
    Sfdrbd   = Sfdrbd + sfdrbd;
    Sgonz    = Sgonz + sgonz;
    Sdzhang  = Sdzhang + sdzhang;
    Smy      = Smy + smy;    
end   
Sfdbd    = Sfdbd/NR;
Sfdrbd   = Sfdrbd/NR;
Sgonz    = Sgonz/NR;
Sdzhang  = Sdzhang/NR;
Smy      = Smy/NR;


%% Results
legendCell{1} = 'FD-BD';
legendCell{2} = 'FD-RBD';
legendCell{3} = 'Gonzalez2018';
legendCell{4} = 'DZhang2019';
legendCell{5} = 'Prop. Alg.';

% filename = 'Figures\MU_OFDM_SR_vs_SNR_data';
% save(filename);

figure
hold on; grid on
plot(SNR_dB,Sfdbd,'k-o','LineWidth',1.5);
plot(SNR_dB,Sfdrbd,'k:o','LineWidth',1.5);
plot(SNR_dB,Sgonz,'g-x','LineWidth',1.5);
plot(SNR_dB,Sdzhang,'b-o','LineWidth',1.5);
plot(SNR_dB,Smy,'r-o','LineWidth',1.5);
legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)','interpreter','latex')
xlabel('SNR (dB)','interpreter','latex')
% filenameFigFile = 'Figures\MU_OFDM_SE_vs_SNR_fig';
% saveas(gcf,filenameFigFile,'pdf');
