%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               CHANNEL DESIGN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to design the mmWave wireless channel matrix.
% Input:
%       Nt - # of antenna elements at transmitter
%       Nr - # of antenna elements at receiver
%       K  - # of users
%       F  - # of subcarrier
%       Ns - # of streams per user
%       Ncl - # of cluster (multipath)
%       Nray - # of rays in each cluster (multipath)
%
% Output: 
%       H  - (K x Nr x Nt) channel matrix
%       At - (K x Ncl*Nray x Nt) channel matrix
%       Ar - (K x Nr x Nt) channel matrix
%       Fopt - (K x Nr x Nt) channel matrix
%       Wopt - (K x Nr x Nt) channel matrix
%
% Description: 
%       The channel consists of the sum of the constribution of Ncl clusters,
%   each cluster consisting of Nray multipath.
%       The AoA and AoD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Guilherme M. Zilli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H,At,Ar,Fopt,Wopt] = channel_realization(Nt,Nr,K,M,Ns,Ncl,Nray,arrayType,Pt)

AoD_Az_ClusterSpread = pi;  % Azimuth cluster mean AoD spread (in rad) 
                            % Uniformly dist. (-AoD_Az_ClusterSpread , AoD_Az_ClusterSpread)
AoD_El_ClusterSpread = pi;  % Elevation cluster mean AoD spread (in rad) 
                            % Uniformly dist. (-AoD_El_ClusterSpread , AoD_El_ClusterSpread)                            
                            
AoA_Az_ClusterSpread = pi;  % Azimuth cluster mean AoD spread (in rad)
                            % Uniformly dist. (-AoA_Az_ClusterSpread , AoA_Az_ClusterSpread)
AoA_El_ClusterSpread = pi;  % Elevation cluster mean AoD spread (in rad)
                            % Uniformly dist. (-AoA_El_ClusterSpread , AoA_El_ClusterSpread)
                            
AoD_Az_VarSpread = 10;     % Azimuth ray angle spread (in deg) - Laplacian dist.
AoD_El_VarSpread = 10;     % Elevation ray angle spread (in deg) - Laplacian dist.

AoA_Az_VarSpread = 10;     % Azimuth ray angle spread (in deg) - Laplacian dist.
AoA_El_VarSpread = 10;     % Elevation ray angle spread (in deg) - Laplacian dist.

    % For ULA array, Elevation Cluster Spread must be set to pi/2 and
    % Elevation Var Spread must be set to zero.

    At = zeros(Nt,Ncl*Nray,K);
    Ar = zeros(Nr,Ncl*Nray,K);
    
%     for k = 1:K
        H = zeros(Nr,Nt,M);
        Fopt = zeros(Nt,Ns,M);
        Wopt = zeros(Nr,Ns,M);
%     end
     
%     for k =1:K
    k=1;
        
        % Angle of Departure
            AoD_Az              = diag(2*(rand(1,Ncl) - 0.5)*AoD_Az_ClusterSpread)*ones(Ncl,Nray);
            AoD_Az_VarSpread_rd = AoD_Az_VarSpread*pi/180;
            AoD_Az              = AoD_Az + laprnd(Ncl,Nray,0,AoD_Az_VarSpread_rd);

            AoD_El              = diag(2*(rand(1,Ncl) - 0.5)*AoD_El_ClusterSpread)*ones(Ncl,Nray);
            AoD_El_VarSpread_rd = AoD_El_VarSpread*pi/180;
            AoD_El              = AoD_El + laprnd(Ncl,Nray,0,AoD_El_VarSpread_rd);

        % Angle of Arrival
            AoA_Az              = diag(2*(rand(1,Ncl) - 0.5)*AoA_Az_ClusterSpread)*ones(Ncl,Nray);
            AoA_Az_VarSpread_rd = AoA_Az_VarSpread*pi/180;
            AoA_Az              = AoA_Az + laprnd(Ncl,Nray,0,AoA_Az_VarSpread_rd);

            AoA_El              = diag(2*(rand(1,Ncl) - 0.5)*AoA_El_ClusterSpread)*ones(Ncl,Nray);
            AoA_El_VarSpread_rd = AoA_El_VarSpread*pi/180;
            AoA_El              = AoA_El + laprnd(Ncl,Nray,0,AoA_El_VarSpread_rd);


        % Normalization factors
            gamma = sqrt((Nt*Nr)/(Ncl*Nray));
            alpha = (1/sqrt(2))*(randn(Ncl,Nray) + 1i*randn(Ncl,Nray));

        % Channel Generation
        for m = 1:M
            for cl = 1:Ncl
                for ray = 1:Nray
                    At(:,(cl-1)*Nray+ray,k) = array_response(arrayType,AoD_Az(cl,ray),AoD_El(cl,ray),Nt);
                    Ar(:,(cl-1)*Nray+ray,k) = array_response(arrayType,AoA_Az(cl,ray),AoA_El(cl,ray),Nr);
                    H(:,:,m) = H(:,:,m) + (gamma*alpha(cl,ray)*reshape(Ar(:,(cl-1)*Nray+ray,k),Nr,1)*reshape(At(:,(cl-1)*Nray+ray,k),Nt,1)')*exp(-1i*2*pi*(cl-1)*(m-1)/M);
                end
            end
        end
        
        % Channel Generation for all subcarriers        
        DELTA = zeros(Nr,Nt);
        DELTA(1:min(Nr,Nt),1:min(Nr,Nt)) = eps*eye(min(Nr,Nt));        
        for m = 1:M
%             if(rank(H{k,m}(:,:))>=Ns)
                [U,S,V] = svd(H(:,:,m) + DELTA); 
                Fopt(:,:,m) = sqrt(Pt) * V(1:Nt,1:Ns)/norm(V(1:Nt,1:Ns),'fro');
                Wopt(:,:,m) = U(1:Nr,1:Ns);
%             else
%                 error('Rank Deficient Channel')
%             end 
        end

%     end
end