function Rsum = MUsumRate(H,F,W,Ns,snr,sigma2_n)
    [Nr,Nt,M,K] = size(H);
    Rsum = 0;
    P = snr*sigma2_n; % just in case sigma2_n =/= 1;
    
    for m = 1:M
        RsumM = 0;
        for k = 1:K
            RsumK = 0;
            for s = 1:Ns
                % INTER-user-interference power
                    interPower = 0;
                    for ki = setdiff(1:K, k)
                        for si = 1:Ns
                            interPower = interPower + abs( W(:,s,m,k)'*H(:,:,m,k)* F(:,si,m,ki) )^2;
                        end
                    end

                % INTRA-user-interference power
                    intraPower = 0;
                    for si = setdiff(1:Ns, s)
                        intraPower = intraPower + abs( W(:,s,m,k)'*H(:,:,m,k)*F(:,si,m,k) )^2;  
                    end

                sinr = ( (snr/(K*Ns))*(abs( W(:,s,m,k)'*H(:,:,m,k)*F(:,s,m,k) )^2) )/( (snr/(K*Ns))*(interPower + intraPower) + sigma2_n*(W(:,s,m,k)'*W(:,s,m,k)) ); 
                RsumK = RsumK + log2(1+sinr);                
            end
            RsumM = RsumM + RsumK;
        end
        Rsum = Rsum + RsumM/M;
    end

end