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
    SNR_dB  = 5;                    % SNR values in dB
    SNR     = 10.^(SNR_dB./10);     % SNR values in decimal

% Scenario Parameters
        K      = 4;            % # of users
        Ns     = 4;                % # or data stream (per user)
        MM     = 2.^[2:9];

    % Transmitter side
        Nt      = 144;              % # of antennas at transmitter
%         Nt_rf   = K*Ns;           % # of RF chains at transmitter
    
    % Receiver side
        Nr      = 16;               % # of antennas at receiver
        Nr_rf   = Ns;               % # of RF chains at receiver
    
    % Channel Parameters
        Ncl     = 5;                % # of clusters
        Nray    = 10;               % # of rays in each cluster
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Pt = K*Ns;
        sigma2 = 0.000000001;
        rho = sigma2*SNR/Pt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Array Type
        arrayType = 'square';
        
        if ( strcmp(arrayType,'square') && ~( (~( mod(Nt,sqrt(Nt)) )) || (~( mod(Nr,sqrt(Nr)) )) ))
            error('The number of antennas must be N = sqrt(N)*sqrt(N) with sqrt(N) integer when using USPA')
        end            
    
    Sfdbd    = zeros(length(MM),1);
    Sfdrbd    = zeros(length(MM),1);
    Sgonz    = zeros(length(MM),1);
    Sdzhang  = zeros(length(MM),1);
    Smy      = zeros(length(MM),1);
    
% Algorithm Simulations
for nreal = 1:NR
    nreal
    
    sfdbd    = zeros(length(MM),1);
    sfdrbd    = zeros(length(MM),1);
    sgonz    = zeros(length(MM),1);
    sdzhang  = zeros(length(MM),1);
    smy      = zeros(length(MM),1);
 
    for mInd = 1:length(MM)
        M = MM(mInd);

        Nt_rf   = K*Ns; 	% # of RF chains at transmitter
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
                
            [ F_fdbd, W_fdbd] = FDBD(H,Ns,Pt);
            
            [ F_fdrbd, W_fdrbd ] = FDRBD(H,Ns,Pt,sigma2,rho);
           
            [ F_Gonz, W_Gonz] = GonzalesComa2018D1( H, Ns, sigma2 ,Pt);

            [ F_DZhang, W_DZhang] = DZhang2019_TensorUnfolding( H, Ns, sigma2 ,Pt);

            [ F_my, W_my] = myMUOFDMBeamformingRCD(H,Ns,sigma2 ,Pt, rho);
                
    
        % Spectral Efficiency
            sfdbd(mInd)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_fdbd,W_fdbd,rho,'total',Pt),W_fdbd,Ns,rho,sigma2);
            sfdrbd(mInd)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_fdrbd,W_fdrbd,rho,'total',Pt),W_fdrbd,Ns,rho,sigma2);
            sgonz(mInd)   = MUspectralEfficiency(H,OFDMWaterFilling(H,F_Gonz,W_Gonz,rho,'total',Pt),W_Gonz,Ns,rho,sigma2);
            sdzhang(mInd) = MUspectralEfficiency(H,OFDMWaterFilling(H,F_DZhang,W_DZhang,rho,'total',Pt),W_DZhang,Ns,rho,sigma2);
            smy(mInd)     = MUspectralEfficiency(H,OFDMWaterFilling(H,F_my,W_my,rho,'total',Pt),W_my,Ns,rho,sigma2);
    end
    
    Sfdbd    = Sfdbd + sfdbd;
    Sfdrbd    = Sfdrbd + sfdrbd;
    Sgonz    = Sgonz + sgonz;
    Sdzhang  = Sdzhang + sdzhang;
    Smy      = Smy + smy;
    
end    
Sfdbd    = Sfdbd/NR;
Sfdrbd    = Sfdrbd/NR;
Sgonz    = Sgonz/NR;
Sdzhang  = Sdzhang/NR;
Smy      = Smy/NR;

%% Results
legendCell{1} = 'FD-BD';
legendCell{2} = 'FD-RBD';
legendCell{3} = 'Gonzalez2018';
legendCell{4} = 'DZhang2019';
legendCell{5} = 'Prop. Alg.';

% filename = strcat('Figures\MU_OFDM_SR_vs_M_data',num2str(SNR_dB));
% save(filename);

figure
hold on; grid on
plot(MM,Sfdbd,'k-o','LineWidth',1.5);
plot(MM,Sfdbd,'k:o','LineWidth',1.5);
plot(MM,Sgonz,'g-x','LineWidth',1.5);
plot(MM,Sdzhang,'b-o','LineWidth',1.5);
plot(MM,Smy,'r-o','LineWidth',1.5);

legend(legendCell,'Location','northwest')
ylabel('Spectral Efficiency (bits/s/Hz)','interpreter','latex')
xlabel('Number of subcarrier - $M$','interpreter','latex')

% filenameFigFile = strcat('Figures\MU_OFDM_SE_vs_M_fig',num2str(SNR_dB));
% saveas(gcf,filenameFigFile,'pdf');