function [F,W] = FDBD(H,Ns,Pt)

    [Nr,Nt,M,K] = size(H);    
%     TOL = 1e-6;
    for m = 1:M
        for j = 1:K
            H_tildeJ = [];
            for k = setdiff(1:K,j)
                H_tildeJ = [H_tildeJ(:,:);   H(:,:,m,k)]; 
            end
            H_tilde(:,:,j) = H_tildeJ;

            [~,St,Vtilde] = svd(H_tilde(:,:,j));
            rt = nnz(St); % Using the nnz function may lead to error since some singular value are extremely small yet nonzero.

%             rt = nnz( max(diag(St),TOL) - TOL);
            if (rt == Nt)
                error('H_tilde is full-rank. Cannot compute null-space')
            end        

            [U,S,V] = svd(H(:,:,m,j)*Vtilde(:,rt+1:end));
            r = nnz(S);

            F(:,:,m,j) = Vtilde(:,rt+1:end)*V(:,1:Ns);
            F(:,:,m,j) = sqrt(Pt/K)* F(:,:,m,j)/norm(F(:,:,m,j),'fro');
            W(:,:,m,j) = U(:,1:Ns);
        end
    end
    
end