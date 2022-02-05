function [F] = OFDMWaterFilling(H,F,W,rho,policy,Pt)
    
    [Nr,Nt,M,K] = size(H);
    [~,Ns,~] = size(F);
    
    switch lower(policy)
        case 'total'
            power = Pt;            
            for m = 1:M
                gain = [];
                for k = 1:K
                    gain = [gain; (abs(diag( W(:,:,m,k)'*H(:,:,m,k)*F(:,:,m,k))))];
%                     gain = [gain; (sqrt(rho)*abs(diag( W(:,:,m,k)'*H(:,:,m,k)*F(:,:,m,k))))];
%                     gain = [gain; ((rho)*abs(diag( W(:,:,m,k)'*H(:,:,m,k)*F(:,:,m,k))).^2)];
                end            
                [ P , ~] = waterfill( 1./gain , power );
                for k = 1:K
                    F(:,:,m,k) = F(:,:,m,k) * diag(sqrt( P(1, (k-1)*Ns + [1:Ns]) ));
                end
            end

        case 'user'
            power = Ns;
            for m = 1:M
                for k = 1:K
                    gain = (rho*abs(diag( W(:,:,m,k)'*H(:,:,m,k)*F(:,:,m,k) )).^2);
                    [ P , ~] = waterfill( 1./gain , power );
                    F(:,:,m,k) = F(:,:,m,k) * diag(sqrt(P));
                end
            end
            
        otherwise
            error('Policy not implemented!')
    end


    function [p, level] = waterfill (x,P)
        % Author: Kenneth Shum, 2010
        % Reference:
        %   T. M. Cover and J. A. Thomas, "Elements of Information Theory", John Wiley & Sons, 2003.
        L=length(x);  % number of channels
        y = sort(x);     % sort the interference power in increasing order
        [a b]= size(y);
        if (a>b)
            y = y' ;  % convert x and y to row vector if necessary
            x = x';
        end
        delta = [0 cumsum((y(2:L)-y(1:(L-1))).*(1:L-1))];
        l = sum(P>=delta); % no. of channel with positive power
        level = y(l)+(P- delta(l))/l;   % water level
        p = (abs(level-x)+level-x)/2;   % the result of pouring water.
    end

end