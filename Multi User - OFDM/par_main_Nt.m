%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   HYBRID BEAMFORMING SIMULATIONS
%   Author: Guilherme Martignago Zilli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Simulation Parameters
    NRlz    = 500;                      % Number of realizations  
    SNR_dB  = 5;                       % SNR values in dB
    SNR     = 10.^(SNR_dB./10);         % SNR values in decimal

% Scenario Parameters
        K       = 4;                    % # of users (K = 1 for SU, K>1 for MU) 
        Ns      = 4;                    % # or data stream (per user)
        M       = 32;

    % Transmitter side
        NT      = [6:16].^2;            % # of antennas at transmitter
        Nt_rf   = K*Ns;                 % # of RF chains at transmitter
    
    % Receiver side
        Nr      = 9;                    % # of antennas at receiver
        Nr_rf   = Ns;                 	% # of RF chains at receiver
    
    % Channel Parameters
        Ncl     = 5;                    % # of clusters
        Nray    = 10;                   % # of rays in each cluster
        
    % Array Type
        arrayType = 'square';
        
    % Transmitted Power
        Pt = K*Ns;
        sigma2 = 0.000000001;
        rho = sigma2*SNR/Pt;
        
    Sfdbd    = zeros(length(NT),1);
    Sfdrbd    = zeros(length(NT),1);
    Sgonz    = zeros(length(NT),1);
    Sdzhang  = zeros(length(NT),1);
    Smy      = zeros(length(NT),1);
 
% Algorithm Simulations
parfor nreal = 1:NRlz
    nreal    
    
    sfdbd    = zeros(length(NT),1);
    sfdrbd    = zeros(length(NT),1);
    sgonz    = zeros(length(NT),1);
    sdzhang  = zeros(length(NT),1);
    smy      = zeros(length(NT),1);
  
    for ntInd = 1:length(NT)
        Nt = NT(ntInd);

        % Checking conditions:
        if ((K*Ns > Nt_rf) || (Nt_rf > Nt))
            error('Check conditions at BS')
        end
        if ((Ns > Nr_rf) || (Nr_rf > Nr))
            error('Check conditions at MS')
        end
        
        % Channel Matrix
        [H,~,~,~,~,~] = channel_realization(Nt,Nr,K,M,Ncl,Nray,arrayType);
        
        % Precoder & Combiner Design
                
            [ F_fdbd, W_fdbd ] = FDBD(H,Ns, Pt); 
            
            [ F_fdrbd, W_fdrbd ] = FDRBD(H,Ns,Pt,sigma2,rho);
           
            [ F_Gonz, W_Gonz ] = GonzalesComa2018D1( H, Ns, sigma2, Pt);

            [ F_DZhang, W_DZhang ] = DZhang2019_TensorUnfolding( H, Ns, sigma2, Pt, rho);

            [ F_my, W_my ] = myMUOFDMBeamformingRCD(H, Ns, sigma2, Pt, rho);
           
        % Spectral Efficiency
            sfdbd(ntInd)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_fdbd,W_fdbd,rho,'total',Pt),W_fdbd,Ns,rho,sigma2);
            sfdrbd(ntInd)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_fdrbd,W_fdrbd,rho,'total',Pt),W_fdrbd,Ns,rho,sigma2);
            sgonz(ntInd)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_Gonz,W_Gonz,rho,'total',Pt),W_Gonz,Ns,rho,sigma2);
            sdzhang(ntInd) = MUspectralEfficiency(H,OFDMWaterFilling(H,F_DZhang,W_DZhang,rho,'total',Pt),W_DZhang,Ns,rho,sigma2);
            smy(ntInd)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);
    end
    
    Sfdbd    = Sfdbd + sfdbd;
    Sfdrbd    = Sfdrbd + sfdrbd;
    Sgonz    = Sgonz + sgonz;
    Sdzhang  = Sdzhang + sdzhang;
    Smy      = Smy + smy;    
end  

Sfdbd    = Sfdbd/NRlz;
Sfdrbd    = Sfdrbd/NRlz;
Sgonz    = Sgonz/NRlz;
Sdzhang  = Sdzhang/NRlz;
Smy      = Smy/NRlz;

%% Results
legendCell{1} = 'FD-BD';
legendCell{2} = 'FD-RBD';
legendCell{3} = 'Gonzalez2018';
legendCell{4} = 'DZhang2019';
legendCell{5} = 'Prop. Alg.';

% filename = strcat('Figures\MU_OFDM_SR_vs_Nt_data',num2str(SNR_dB));
% save(filename);

figure
hold on; grid on
plot(NT,Sfdbd,'k-o','LineWidth',1.5);
plot(NT,Sfdrbd,'k:o','LineWidth',1.5);
plot(NT,Sgonz,'g-x','LineWidth',1.5);
plot(NT,Sdzhang,'b-o','LineWidth',1.5);
plot(NT,Smy,'r-o','LineWidth',1.5);

legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)','interpreter','latex')
xlabel('Number of transmit antennas - $N_t$','interpreter','latex')

% filenameFigFile = strcat('Figures\MU_OFDM_SE_vs_Nt_fig',num2str(SNR_dB));
% saveas(gcf,filenameFigFile,'pdf');
