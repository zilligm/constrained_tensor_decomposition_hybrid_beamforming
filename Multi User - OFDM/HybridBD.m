function [F,W] = HybridBD(H,Ns)

    [Nr,Nt,K] = size(H);    
    HS = [];
   
    for j = 1:K
        jNot = 1:K;
        jNot(j) = [];
        H_tildeJ = [];
        for k = jNot
            H_tildeJ = [H_tildeJ(:,:); ...
                              reshape(H(:,:,k),Nr,Nt)]; 
        end
        H_tilde(:,:,j) = H_tildeJ;
        
        [~,St,Vtilde] = svd(squeeze(H_tilde(:,:,j)));
        rt = nnz(St);
        if (rt == Nt)
            error('H_tilde is full-rank. Cannot compute null-space')
        end
        
        
        [U,S,V] = svd(reshape(H(:,:,j),Nr,Nt)*Vtilde(:,rt+1:end));   % Water-filling not implemented.If required, I must use the singlular values here.
        r = nnz(S);
               
        F(:,:,j) = Vtilde(:,rt+1:end)*V(:,1:Ns);
        W(:,:,j) = U(:,1:Ns);
    end
    
end