function [ F, W] = DZhang2019_TensorUnfolding( H, Ns, sigma2, Pt, rho)
    if length(size(H)) == 4
        [Nr,Nt,M,K] = size(H);
    else
        error('Invalid channel dimensions');
    end      

    
    % Compute analog precoder and combiner
        WRF = zeros(Nr,Ns,K);
        FRF = zeros(Nt,K*Ns);        
   
        H1 = zeros(Nr,Nt*M);
        
        for k = 1:K 
            H1 = zeros(Nr,Nt*M);
            for m = 1:M
                H1(:, (m-1)*Nt + 1: m*Nt) = H(:,:,m,k);
            end
            
            [U,~,~] = svds(H1,Ns);
            WRF(:,:,k) = (1/sqrt(Nr))*exp(1i*angle(U(:,1:Ns)));
        end
        
        for k = 1:K 
            H1 = zeros(M*Ns,Nt);
            for m = 1:M
                H1((m-1)*Ns + 1: m*Ns, : ) = WRF(:,:,k)' * H(:,:,m,k);
            end
            
            [~,~,V] = svds(H1,Ns);
            FRF(:, (k-1)*Ns + [1:Ns] ) = (1/sqrt(Nt))*exp(1i*angle(V(:,1:Ns)));
        end


    
    % Compute the baseband equivalent channel
        for m = 1:M
            for k = 1:K
                Heq(:,:,m,k) = WRF(:,:,k)'*H(:,:,m,k)*FRF;
            end
        end
        
    % Compute digital BF   
        for m = 1:M
            for k = 1:K            
                H_tildeJ = [];
                for j = setdiff(1:K,k)
                    H_tildeJ = [H_tildeJ(:,:); Heq(:,:,m,j)]; 
                end   
                
                [V,~] = eig( Heq(:,:,m,k)'*Heq(:,:,m,k) , (K*(Ns^2)*sigma2/Pt)*eye(K*Ns) + H_tildeJ'*H_tildeJ );
%                 [V,~] = eig( Heq(:,:,m,k)'*Heq(:,:,m,k) , (K*Ns^2*sigma2/(rho*Pt))*eye(K*Ns) + H_tildeJ'*H_tildeJ );
%                 if K == 7       %% This fixes a Matlab bug that I couldn't figure out 
%                     V = V(:,1:Ns);
%                 else
%                     V = V(:,end:-1:end-Ns+1);
%                 end
                
                [~,Ind] = sort(sum(abs(Heq(:,:,m,k)*V),1),'descend');
                V = V(:,Ind(1:Ns));
                
                Fbb(:,:,m,k) = sqrt(Pt/K)* V /norm(FRF*V,'fro');
                Wbb(:,:,m,k) = pinv(Heq(:,:,m,k)*Fbb(:,:,m,k)*Fbb(:,:,m,k)'*Heq(:,:,m,k)' + (K*Ns*sigma2/(Pt))*WRF(:,:,k)'*WRF(:,:,k) )*Heq(:,:,m,k)*Fbb(:,:,m,k);
%                 Wbb(:,:,m,k) = pinv(Heq(:,:,m,k)*Fbb(:,:,m,k)*Fbb(:,:,m,k)'*Heq(:,:,m,k)' + (K*Ns*sigma2/(rho*Pt))*WRF(:,:,k)'*WRF(:,:,k) )*Heq(:,:,m,k)*Fbb(:,:,m,k);
                Wbb(:,:,m,k) = sqrt(Ns)* Wbb(:,:,m,k) /norm( WRF(:,:,k)*Wbb(:,:,m,k),'fro');
                
            end 
        
            for k = 1:K
%                 Fbb(:,:,m,j) = b*Fbb(:,:,m,j);
                F(:,:,m,k) = FRF*Fbb(:,:,m,k);
                W(:,:,m,k) = WRF(:,:,k)*Wbb(:,:,m,k);
            end
        end
end