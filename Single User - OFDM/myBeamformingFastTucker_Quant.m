function [ F, W, DATA] = myBeamformingFastTucker2_Quant( H, Ns, Nb, Pt)
    [Nr,Nt,M] = size(H); 

    W = zeros(Nr,Ns,M);
    F = zeros(Nt,Ns,M);
    
    WRF = zeros(Nr,Ns);
    FRF = zeros(Nt,Ns);
        
    epsilon = 0.001;
    maxIte  = 20;
    
    DATA = zeros(Ns,1);
    
    res = 2*pi/(2^Nb);
    quantAngle = 0:res:2*pi-res;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Hres = H;
    
    for s = 1:Ns
        p = complex(randn(Nr,1),randn(Nr,1));   % Randomly initialization
        p = (1/sqrt(Nr))*p./abs(p);
        q = complex(randn(Nt,1),randn(Nt,1));   % Randomly initialization
        q = (1/sqrt(Nt))*q./abs(q);
            
        count = 0;
%         delta = [2*epsilon 0];
        delta = [2*epsilon 0 zeros(1,maxIte)];
        
%         while ((count < maxIte)&&(abs(delta(end)-delta(end-1))> epsilon))
        while ((count < maxIte)&&(abs(delta(count+2)-delta(count+1))> epsilon))
            count = count + 1;
            
%             HresP = zeros(Nr,Nr);
%             for m = 1:M
%                 a = Hres(:,:,m)*q;
%                 HresP = HresP + a*a';
%             end
%             p = HresP*p;
%             p = (1/sqrt(Nr))*p./abs(p); 
%             
%             HresQ = zeros(Nt,Nt);
%             for m = 1:M
%                 b = Hres(:,:,m)'*p;
%                 HresQ = HresQ + b*b';
%             end
%             q = HresQ*q;
%             q = (1/sqrt(Nt))*q./abs(q);

            pp = zeros(Nr,1);
            for m = 1:M
                a = Hres(:,:,m)*q;
                b = a'*p;
                pp = pp + a*b;
            end
            p = pp;
%             p = (1/sqrt(Nr))*p./abs(p); 
            p = PSquant( (1/sqrt(Nr))*p./abs(p) , quantAngle);

            
            qq = zeros(Nt,1);
            HresX = zeros(Nt,Nt);
            for m = 1:M
                a = Hres(:,:,m)'*p;
                HresX = HresX + a*a';
                b = a'*q;
                qq = qq + a*b;
            end
            q = qq;
%             q = (1/sqrt(Nt))*q./abs(q);
            q = PSquant( (1/sqrt(Nt))*q./abs(q) , quantAngle);

%             delta = [delta abs(q'*HresQ*q)/M];
            delta(count+2) = abs(q'*HresX*q)/M;
        end

        WRF(:,s) = p;
        FRF(:,s) = q;
        
%         PrL = (eye(Nr) - p*p');
%         PrR = (eye(Nt) - q*q');
% 
%         for m = 1:M 
%             Hres(:,:,m) = PrL * Hres(:,:,m) * PrR;
%         end   

        for m = 1:M 
            a = p'*Hres(:,:,m);
            Hres(:,:,m) = Hres(:,:,m) - p*a;
            b = Hres(:,:,m)*q; 
            Hres(:,:,m) = Hres(:,:,m) - b*q';
        end 
        
        DATA(s,1) = count;
    end    
%     FRF = Q;
%     WRF = P;
    for m = 1:M
        Heff = WRF'*H(:,:,m)*FRF;
        [U,~,V] = svd(Heff);
%         FBB_my = ;
%         FBB_my = sqrt(Ns) * V(:,1:Ns) / norm(FRF_my*V(:,1:Ns),'fro');
        F(:,:,m) = sqrt(Pt) * FRF * V(:,1:Ns) / norm(FRF*V(:,1:Ns),'fro');
%         W = U(:,1:Ns);
        W(:,:,m) = WRF*U(:,1:Ns);
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