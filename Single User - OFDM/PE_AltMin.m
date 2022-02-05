function [FRF, FBB, count] = PE_AltMin(Fopt, NRF);
      [Nt,Ns,M] = size(Fopt);
%     [Nr,Nt,M] = size(H);
    
%     if (K ~= 1)
%         error('This algorithm was designed solely for the Single User case.');
%     end

    mynorm = [];
    FRF = exp( sqrt(-1) * unifrnd (0,2*pi,Nt,NRF) );
    %
    count = 0;
    %
    while (isempty(mynorm) || abs( mynorm(end) - mynorm(end-1) ) > 1e-2)
        AUX = zeros(Nt,Ns);
        auxConv = 0;
        for m = 1:M
            [U,S,V] = svd(Fopt(:,:,m)'*FRF);
            FBB(:,:,m) = V(:,[1:Ns])*U';
            
            AUX = AUX + Fopt(:,:,m)*FBB(:,:,m)';
            
            auxConv = auxConv + norm(Fopt(:,:,m) * FBB(:,:,m)' - FRF,'fro')^2;
        end
        mynorm = [mynorm, auxConv/M];
        
        FRF = (1/sqrt(Nt))*exp(1i * angle(AUX));
        
        auxConv = 0;
        for m = 1:M
            auxConv = auxConv + norm(Fopt(:,:,m) * FBB(:,:,m)' - FRF,'fro')^2;
        end
        mynorm = [mynorm, auxConv/M];
        
        count = count + 1;

    end
    

end
