function [F,W] = THTsai2019( H, Nt_rf, Nr_rf, Ns, Pt)
    [Nr,Nt,M] = size(H);

    % Analog Combiner Design
        HS1 = zeros(Nr,Nr);
        for m = 1:M
            HS1 = HS1 + H(:,:,m)*H(:,:,m)';
        end
        [US1,SS1,VS1] = svd(HS1);
        WRF = US1(:,1:Nr_rf)./abs(US1(:,1:Nr_rf));

    % Analog Precoder Design
        HS2 = zeros(Nt,Nt);
        for m = 1:M
            HS2 = HS2 + H(:,:,m)'*H(:,:,m);
        end
        [US2,SS2,VS2] = svd(HS2);
        FRF = US2(:,1:Nt_rf)./abs(US2(:,1:Nt_rf));
        
    %%% Digital Precoder and Combiner Design
        for m = 1:M
            Heff = WRF'*H(:,:,m)*FRF;
            [Ue,Se,Ve] = svd(Heff);
            
            FBB(:,:,m) = Ve(:,1:Ns);
            FBB(:,:,m) = sqrt(Pt) * FBB(:,:,m) / norm(FRF * FBB(:,:,m),'fro');
            F(:,:,m) = FRF*FBB(:,:,m); % Effective Precoder
            
            WBB(:,:,m) = Ue(:,1:Ns);            
            W(:,:,m) = WRF*WBB(:,:,m); % Effective Combiner
        end

end