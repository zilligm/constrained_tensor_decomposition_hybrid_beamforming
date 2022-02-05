function [ F, W, DATA] = Chiang2018( H, Ns, CDBK, S, sigma2, Pt )
    [Nr,Nt,M] = size(H); 
    [~,Nc] = size(CDBK);
    
%     S = S; % This is the M in the paper (called S here since M is the number of subcarriers)
    Ic = nchoosek(S,Ns);
    
    Y = zeros(Nc,Nc,M);
    YY = zeros(Nc,Nc);
    for m = 1:M
        Y(:,:,m) = CDBK'*H(:,:,m)*CDBK + ((1/sqrt(2))*sigma2)*(randn(Nc,Nc)+1i*randn(Nc,Nc));
        YY = YY + abs(Y(:,:,m)).^2;
    end
    
    Icom = [];
    Ipre = [];
    % Selecting candidates
    for s = 1:S
    	[com,pre] = find(YY == max(max(YY)));
        Ipre = [Ipre pre];
        Icom = [Icom com];
        YY(:,pre) = 0;
        YY(com,:) = 0;
    end
    
    %%%%%%%%%%%%%%%%%%% 
    % Needs to be fixed to include the case where S>Ns
    % Here, I'd solve (26)
    if S == Ns
        FP = CDBK(:,Ipre);
        WP = CDBK(:,Icom);
    else
        If = nchoosek(Ipre,Ns);
        Iw = nchoosek(Icom,Ns);
        
        [Ni,~] = size(If);
        
        HFro = zeros(Ni,Ni);
        
        for m = 1:M
            
            for ni1 = 1:Ni
                for ni2 = 1:Ni
%                     HFro(ni1,ni2) = norm( Y(Iw(ni1,:),If(ni2,:),m) ,'fro')^2;
                    HFro(ni1,ni2) = det( Y(Iw(ni1,:),If(ni2,:),m))^2;
                end
            end
            
        end
        
        [iW,iF] = find(HFro == max(max(HFro)));
        Icom = Iw(iW,:);
        Ipre = If(iF,:);
        FP = CDBK(:,Ipre);
        WP = CDBK(:,Icom);
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%
    for m = 1:M
%         HE_hat = WP' * H(:,:,m) * FP;
        HE_hat = Y(Icom,Ipre,m);
        He(:,:,m) = (WP'*WP)^(-0.5) * ( HE_hat ) * (FP'*FP)^(-0.5);
%         He(:,:,m) = pinv( sqrt(WP'*WP) ) * ( HE_hat ) * pinv( sqrt(FP'*FP) );
        [Ue,~,Ve] = svd(He(:,:,m));
        
        FB = (FP'*FP)^(-0.5)*Ve(:,1:Ns);
        FB = sqrt(Pt) * FB / norm(FP*FB,'fro');
        WB = (WP'*WP)^(-0.5)*Ue(:,1:Ns);
        
        F(:,:,m) = FP*FB;
        W(:,:,m) = WP*WB;
    end

end