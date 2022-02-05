function [ F, W] = GonzalesComa2018D1( H, Ns, sigma2, Pt, rho)
    if length(size(H)) == 4
        [Nr,Nt,M,K] = size(H);
    else
        error('Invalid channel dimensions');
    end      

    % Compute analog precoder and combiner
        WRF = zeros(Nr,Ns,K);
        FRF = zeros(Nt,K*Ns);   
        
        WFD = zeros(Nr,Ns,M,K);
        FFD = zeros(Nt,Ns,M,K);
            
        % Fully digital DL combiner (UL precoder)
            for m = 1:M
                for k = 1:K
                	[~,~,V] = svd(H(:,:,m,k)');
                    WFD(:,:,m,k) = V(:,1:Ns);
                end                
            end 
            
        % Fully digital DL precoder (UL combiner)
                % MMSE
                for m = 1:M
                    for k = 1:K
%                         R = (sigma2)*eye(Nt);
                        R = (sigma2)*eye(Nt);
                        for j = setdiff(1:K,k)
                        R = R + H(:,:,m,j)'*WFD(:,:,m,j)*WFD(:,:,m,j)'*H(:,:,m,j);
                        end
                        FFD(:,:,m,k) = inv(R)*H(:,:,m,k)'*WFD(:,:,m,k);
                    end
                end
                    
                    
                % MRC
%                 for m = 1:M
%                     for k = 1:K
%                         FFD(:,:,m,k) = H(:,:,m,k)'*WFD(:,:,m,k);
%                     end
%                 end                    
                    

                % BD
%                     for m = 1:M
%                         for k = 1:K
%                             H_tildeJ = [];
%                             for j = setdiff(1:K,k)
%                                 H_tildeJ = [H_tildeJ;   H(:,:,m,j)]; 
%                             end
%                             [~,St,Vtilde] = svd(H_tildeJ);
%                             rt = nnz(St); 
%                             if (rt == Nt)
%                                 error('H_tilde is full-rank. Cannot compute null-space')
%                             end        
% 
%                             [U,S,V] = svd(H(:,:,m,k)*Vtilde(:,rt+1:end));
%                             r = nnz(S);
% 
%                             FFD(:,:,m,k) = Vtilde(:,rt+1:end)*V(:,1:Ns);
%                             FFD(:,:,m,k) = sqrt(Ns)* FFD(:,:,m,k)/norm(FFD(:,:,m,k),'fro');
%                             WFD(:,:,m,k) = U(:,1:Ns);
%                         end
%                     end
                
            
        % DL precoder factorization
            [Fbb,FRF] = ProjGradFactorizationF(FFD,Pt);
            
            for k = 1:K                
                [Wbb(:,:,:,k),WRF(:,:,k)] = ProjGradFactorizationF(WFD(:,:,:,k),Pt);
            end
         

        for m = 1:M
            for k = 1:K
                F(:,:,m,k) = sqrt(Pt/K)*FRF*Fbb(:,:,m,k)/norm(FRF*Fbb(:,:,m,k),'fro');
                W(:,:,m,k) = WRF(:,:,k)*Wbb(:,:,m,k);
            end                
        end

            
end

function [ABB,ARF] = ProjGradFactorizationF(AFD,Pt)
    if ndims(AFD) == 4
        [N,Ns,M,K] = size(AFD);
    elseif ndims(AFD) == 3
        [N,Ns,M] = size(AFD);
        K = 1;
    end
    
    ARF = randn(N,K*Ns) + 1i*randn(N,K*Ns);
    ARF = (1/sqrt(N))*ARF./abs(ARF);
    s = 1;
    c = 0;
    maxIte = 100;
    d = 0;
    for m = 1:M
        for k = 1:K
            A(:,:,m,k) = sqrtm(ARF'*ARF)*ARF'*AFD(:,:,m,k);
            d = d + norm(AFD(:,:,m,k),'fro')^2 - 2*sqrt(Pt/Ns)*norm(A(:,:,m,k)) + Pt/Ns;
        end
    end
    
    delta = [Inf d];
    while (abs(delta(end) - delta(end-1))>0.001)&&(c < maxIte)
        c = c+1;
        deld = zeros(N,K*Ns);
        for m = 1:M
            for k = 1:K
                deld = deld + ((Pt/Ns)/norm(A(:,:,m,k)))*(ARF*pinv(ARF'*ARF)*ARF' - eye(N))*AFD(:,:,m,k)*AFD(:,:,m,k)'*ARF*pinv(ARF'*ARF);
            end
        end
        
        ARF = ARF - s*deld;
        ARF = (1/sqrt(N))*ARF./abs(ARF);
        s = s/1.2;
        
        d = 0;
        for m = 1:M
            for k = 1:K
                A(:,:,m,k) = sqrtm(ARF'*ARF)*ARF'*AFD(:,:,m,k);
                d = d + norm(AFD(:,:,m,k),'fro')^2 - 2*sqrt(Pt/Ns)*norm(A(:,:,m,k)) + Pt/Ns;
            end
        end
        delta = [delta d];
    end
        
    for m = 1:M
        for k = 1:K
            ABB(:,:,m,k) = pinv(ARF)*AFD(:,:,m,k);
            ABB(:,:,m,k) = sqrt(Ns)*ABB(:,:,m,k) / norm(ARF*ABB(:,:,m,k),'fro');  
        end
    end
    
end