function [ F, W] = myMUOFDMBeamformingRCD_Quant( H, Ns, sigma2, Pt, rho, Nb )
    if length(size(H)) == 4
        [Nr,Nt,M,K] = size(H);
    else
        error('Invalid channel dimensions');
    end      

    res = 2*pi/(2^Nb);
    quantAngle = 0:res:2*pi-res;
    
    % Compute analog precoder and combiner
        WRF = zeros(Nr,Ns,K);
        FRF = zeros(Nt,K*Ns);        
        epsilon = 0.001;
        maxIte  = 20;    
%         DATA = zeros(Ns,1);
        Hres = H;    
        for s = 1:Ns
            for k = 1:K 
                w = (1/sqrt(Nr))*exp(1i*angle(2*pi*rand(Nr,1))); % Randomly initialization
                f = (1/sqrt(Nt))*exp(1i*angle(2*pi*rand(Nt,1))); % Randomly initialization
                delta = 1;
                delta = [delta 0];        
                count = 0;
                while ((count < maxIte)&&(abs(delta(end)-delta(end-1))> epsilon))
                    count = count + 1;
                    
                    HresP = zeros(Nr,Nr);
                    for m = 1:M
                        HresP = HresP + Hres(:,:,m,k)*f*f'*Hres(:,:,m,k)';
                    end
%                     w = (1/sqrt(Nr))*exp(1i*angle(HresP*w));
                    w = PSquant( (1/sqrt(Nr))*exp(1i*angle(HresP*w)) , quantAngle);

                    HresQ = zeros(Nt,Nt);
                    for m = 1:M
                        HresQ = HresQ + Hres(:,:,m,k)'*w*w'*Hres(:,:,m,k);
                    end
%                     f = (1/sqrt(Nt))*exp(1i*angle(HresQ*f));   
                    f = PSquant( (1/sqrt(Nt))*exp(1i*angle(HresQ*f)) , quantAngle);

                    delta = [delta abs(f'*HresQ*f)];
                                
                end
                
                WRF(:,s,k) = w;
                FRF(:, (k-1)*Ns + s )  = f;
                
                for m = 1:M
                    Hres(:,:,m,k) = (eye(Nr) - w*w') * Hres(:,:,m,k);
%                     Hres(:,:,m,k) = (eye(Nr) - WRF(:,:,k)*pinv(WRF(:,:,k)'*WRF(:,:,k))*WRF(:,:,k)') * Hres(:,:,m,k);
                    for kc = 1:K
                        Hres(:,:,m,kc) = Hres(:,:,m,kc) * (eye(Nt) - f*f');
%                         Hres(:,:,m,kc) = Hres(:,:,m,kc) * (eye(Nt) - FRF*pinv(FRF'*FRF)*FRF');
                    end
                end

                
            end
        end
    
    % Compute the baseband equivalent channel
        for m = 1:M
            for k = 1:K
                Heq(:,:,m,k) = WRF(:,:,k)'*H(:,:,m,k)*FRF;
            end
        end
        
    % Compute digital BF using REGULARIZED CHANNEL DIAGONALIZATION  
        for m = 1:M
            for j = 1:K            
                H_tildeJ = [];
                for k = setdiff(1:K,j)
                    H_tildeJ = [H_tildeJ(:,:); Heq(:,:,m,k)]; 
                end         
                Bk = pinv(rho*H_tildeJ'*H_tildeJ + ((K*Ns*sigma2)/Pt)*eye(K*Ns));
                [U,~,V] = svd(Heq(:,:,m,j)*Bk);

                Fbb(:,:,m,j) = Bk*V(:,1:Ns);            
%                 Fbb(:,:,j) = sqrt(Ns)* Fbb(:,:,j)/norm(Frf*Fbb(:,:,j),'fro');  % Per user power
    %             F(:,:,j) = Frf*Fbb(:,:,j);                                     % Per user power

                Wbb(:,:,m,j) = U(:,1:Ns);
            end 
        
%             B = [];
%             for k = 1:K
%                 B = [B Fbb(:,:,m,k)];
%             end
%             b = sqrt(K*Ns)/norm(FRF*B,'fro');

            for j = 1:K
%                 Fbb(:,:,m,j) = b*Fbb(:,:,m,j);
                Fbb(:,:,m,j) = sqrt(Ns)*Fbb(:,:,m,j)/norm(FRF*Fbb(:,:,m,j),'fro');
                F(:,:,m,j) = FRF*Fbb(:,:,m,j);
                W(:,:,m,j) = WRF(:,:,j)*Wbb(:,:,m,j);
            end
        end
end


function v = PSquant(v,qVal)
    N = length(v);
    for c = 1:N
        ang = angle(v(c));
        if ang < 0
            ang = ang + 2*pi;
        end            
        [~,b] = min(abs(qVal-ang));
        v(c) = (1/sqrt(N))*exp(1i*qVal(b));
    end
end