function R = SUspectralEfficiency(H,F,W,Ns,rho,sigma2)
    [Nr,Nt,M] = size(H);
    
    R = 0;              % Spectral Efficiency Initialization
    
    for m = 1:M
        Rm = 0;         % Spectral Efficiency of the m-th subcarrier;

        Hi = H(:,:,m);
        Fi = F(:,:,m);
        Wi = W(:,:,m);
        
        Rn = Wi'*Wi + eps*eye(Ns);
        R = R + (1/M) * real( log2( det( eye(Ns) + (rho/sigma2)*pinv(Rn)*(Wi'*(Hi*(Fi*Fi')*Hi')*Wi) ) ));
    end
end