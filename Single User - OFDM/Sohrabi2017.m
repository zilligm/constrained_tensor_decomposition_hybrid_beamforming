function [F,W, nite] = Sohrabi2017(H,Ns, sigma2,Pt)

    [Nr,Nt,M] = size(H);
    
%     if (K ~= 1)
%         error('This algorithm was designed solely for the Single User case.');
%     end

    epsilon = 0.01;
    %%% Finding the analog precoder using Algorithm 1 in [15] (Sohrabi2016)
        % This algorithm requires F1, gamma^2, and sigma^2
        % According to [15], gamma^2 = P/(Nt * Ns);
        
        F1 = zeros(Nt,Nt);
        for m = 1:M
            F1 = F1 + H(:,:,m)'*H(:,:,m);
        end
        F1 = F1/M;
        
        FRF = ones(Nt,Ns);

        count = 0;
        maxIte = 100;
        testConv = [0 norm(FRF)];     
        while ( (count <= maxIte) && (testConv(end) >= epsilon) )
            count = count + 1;
            Faux = FRF;
            
            for s = 1:Ns
                FRFbar = FRF(:,setdiff(1:Ns, s));
                Cs = eye(Ns-1) + (Pt/(Nt*Ns*sigma2))*FRFbar'*F1*FRFbar;
                Gs = (Pt/(Nt*Ns*sigma2))*F1 - ((Pt/(Nt*Ns*sigma2))^2)*F1*FRFbar*pinv(Cs)*FRFbar'*F1;
                for n = 1:Nt
                    eta = Gs(n,:)*FRF(:,s) - Gs(n,n)*FRF(n,s);
                    if eta == 0
                        FRF(n,s) = (1/sqrt(Nt))*1;
                    else
                        FRF(n,s) = (1/sqrt(Nt))*eta/abs(eta);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%
%                     if eta == 0
%                         FRF(n,s) = 1;
%                     else
%                         FRF(n,s) = eta/abs(eta);
%                     end
                    %%%%%%%%%%%%%%%%%%%%%%
                end
            end
            testConv = [testConv norm(Faux - FRF)];
        end
        
        nite(1,1) = count; 
     
    %%% Finding the digital precoder
        Q = FRF'*FRF;
        QinvSqrt = pinv(Q^(1/2));
        for m = 1:M
            Heff = H(:,:,m)*FRF;
            [~,~,Ve] = svd(Heff);
            Le =  sqrt(Pt/Ns)*eye(Ns);% This matrix corresponds to the power allocation over the subcarriers; For now, it's assumed equal power allocation.
            FBB(:,:,m) = QinvSqrt*Ve(:,1:Ns)*Le;
            FBB(:,:,m) = sqrt(Pt) * FBB(:,:,m) / norm(FRF * FBB(:,:,m),'fro');
            F(:,:,m) = FRF*FBB(:,:,m); % Effective Precoder
        end
    
    %%% Finding the analog combiner using Algorithm 1 in [15] (Sohrabi2016)
        F2 = zeros(Nr,Nr);
        for m = 1:M
            F2 = F2 + H(:,:,m)*F(:,:,m)*(H(:,:,m)*F(:,:,m))';
        end
        F2 = F2/M;
        
        WRF = ones(Nr,Ns);

        count = 0;
        maxIte = 100;
        testConv = [0 norm(WRF)];
        while ( (count <= maxIte) && (testConv(end) >= epsilon) )
            count = count + 1;
            Waux = WRF;
            
            for s = 1:Ns
                WRFbar = WRF(:,setdiff(1:Ns, s));
                Cs = eye(Ns-1) + (1/(Nr*sigma2))*WRFbar'*F2*WRFbar;
                Gs = (1/(Nr*sigma2))*F2 - ((1/(Nr*sigma2))^2)*F2*WRFbar*pinv(Cs)*WRFbar'*F2;
                for n = 1:Nr
                    eta = Gs(n,:)*WRF(:,s) - Gs(n,n)*WRF(n,s);
                    if eta == 0
                        WRF(n,s) = (1/sqrt(Nr))*1;
                    else
                        WRF(n,s) = (1/sqrt(Nr))*eta/abs(eta);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%
%                     if eta == 0
%                         WRF(n,s) = 1;
%                     else
%                         WRF(n,s) = eta/abs(eta);
%                     end
                    %%%%%%%%%%%%%%%%%%%%
                end
            end
            testConv = [testConv norm(Waux - WRF)];
        end
        nite(2,1) = count;
        
    %%% Finding the digital combiner
%         Q = FRF'*FRF;
%         QinvSqrt = pinv(Q^(1/2));
        for m = 1:M
            AUX = WRF'*H(:,:,m)*F(:,:,m);
            J =   (AUX*AUX') + sigma2*WRF'*WRF;
            WBB(:,:,m) = pinv(J)*AUX;
            
%             WBB(:,:,m)  = sqrt(Pt) * WBB(:,:,m)  / norm(WRF * WBB(:,:,m),'fro');
            W(:,:,m) = WRF*WBB(:,:,m); % Effective Combiner
        end

end