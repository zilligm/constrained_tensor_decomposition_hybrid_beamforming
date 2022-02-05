function [ F, W, DATA] = myBeamformingFastTucker( H, Ns, Pt )
    [Nr,Nt,M] = size(H); 

    W = zeros(Nr,Ns,M);
    F = zeros(Nt,Ns,M);
    
    WRF = zeros(Nr,Ns);
    FRF = zeros(Nt,Ns);
        
    epsilon = 0.001;
    maxIte  = 30;
    
    DATA = zeros(Ns,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hres = H;
    
    for s = 1:Ns
        p = complex(randn(Nr,1),randn(Nr,1));   % Randomly initialization
        p = (1/sqrt(Nr))*p./abs(p);
        q = complex(randn(Nt,1),randn(Nt,1));   % Randomly initialization
        q = (1/sqrt(Nt))*q./abs(q);
            
        count = 0;
        delta = [2*epsilon 0 zeros(1,maxIte)];
        
        while ((count < maxIte)&&(abs(delta(count+2)-delta(count+1))> epsilon))
            count = count + 1;
            
            % For higher efficiency, expression in (18) and (19) in the
            % original paper are broken down as a summation of smaller
            % matrix products (which is possible thanks to the block
            % diagonal structure of kron(I_M, v*v')

            pp = zeros(Nr,1);
            for m = 1:M
                a = Hres(:,:,m)*q;
                b = a'*p;
                pp = pp + a*b;
            end
            p = pp;
            p = (1/sqrt(Nr))*p./abs(p); 
            
            qq = zeros(Nt,1);
            HresX = zeros(Nt,Nt);
            for m = 1:M
                a = Hres(:,:,m)'*p;
                HresX = HresX + a*a';
                b = a'*q;
                qq = qq + a*b;
            end
            q = qq;
            q = (1/sqrt(Nt))*q./abs(q);

%             delta = [delta abs(q'*HresQ*q)/M];
            delta(count+2) = abs(q'*HresX*q)/M;
        end

        WRF(:,s) = p;
        FRF(:,s) = q;

        for m = 1:M 
            a = p'*Hres(:,:,m);
            Hres(:,:,m) = Hres(:,:,m) - p*a;
            b = Hres(:,:,m)*q; 
            Hres(:,:,m) = Hres(:,:,m) - b*q';
        end 
        
        DATA(s,1) = count;
    end    

    for m = 1:M
        Heff = WRF'*H(:,:,m)*FRF;
        [U,~,V] = svd(Heff);
        F(:,:,m) = sqrt(Pt) * FRF * V(:,1:Ns) / norm(FRF*V(:,1:Ns),'fro');
        W(:,:,m) = WRF*U(:,1:Ns);
    end
end