function [ F, W] = JZhang2016( H, Ns, Pt )
    [Nr,Nt,M] = size(H); 

    W = zeros(Nr,Ns,M);
    F = zeros(Nt,Ns,M);
    
    WRF = zeros(Nr,Ns);
    FRF = zeros(Nt,Ns);
        
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U1,~,~] = svds(matrix(H,1),Ns); U1 = U1(:,1:Ns);
    [U2,~,~] = svds(matrix(H,2),Ns); U2 = U2(:,1:Ns);

    WRF = (1/sqrt(Nr))*U1./abs(U1);
    FRF = (1/sqrt(Nt))*conj(U2)./abs(U2);

    for m = 1:M
        Heff = WRF'*H(:,:,m)*FRF;
        [U,~,V] = svd(Heff);
        F(:,:,m) = sqrt(Pt) * FRF * V(:,1:Ns) / norm(FRF*V(:,1:Ns),'fro');
        W(:,:,m) = WRF*U(:,1:Ns);
    end
end



function Y = matrix(X,n)
    [A,B,C] = size(X);
    Y = [];
    switch n
        case 1
            for c = 1:C
                Y = [Y X(:,:,c)];
            end
        case 2
            for c = 1:C
                Y = [Y transpose(X(:,:,c))];
            end
        case 3
            for c = 1:C
                Y = [Y; transpose(vec(X(:,:,c)))];
            end
        otherwise
            disp('Invalid mode')
    end
end

function x = Tnorm(X)
    sizeX = size(X);
    x = 0;
    for i1 = 1:sizeX(1)
        for i2 = 1:sizeX(2)
            for i3 = 1:sizeX(3)
                x = x + abs(X(i1,i2,i3)'*X(i1,i2,i3));
            end
        end
    end
    x = sqrt(x);
end

function Y = nprod(X,A,n)
    sizeX = size(X);
    if length(size(A)) == 1
        J = length(A);
        In = 1;
    else
        [J,In] = size(A);
    end
    if (In == sizeX(n))
        switch n
            case 1
                Y = zeros(J,sizeX(2),sizeX(3));
                for i1 = 1:J
                    for i2 = 1:sizeX(2)
                        for i3 = 1:sizeX(3)
                            for in = 1:In
                                Y(i1,i2,i3) = Y(i1,i2,i3) + X(in,i2,i3)*A(i1,in);
                            end
                        end
                    end
                end
                
            case 2
                Y = zeros(sizeX(1),J,sizeX(3));
                for i1 = 1:sizeX(1)
                    for i2 = 1:J
                        for i3 = 1:sizeX(3)
                            for in = 1:In
                                Y(i1,i2,i3) = Y(i1,i2,i3) + X(i1,in,i3)*A(i2,in);
                            end
                        end
                    end
                end
                
            case 3
                Y = zeros(sizeX(1),sizeX(2),J);
                for i1 = 1:sizeX(1)
                    for i2 = 1:sizeX(2)
                        for i3 = 1:J
                            for in = 1:In
                                Y(i1,i2,i3) = Y(i1,i2,i3) + X(i1,i2,in)*A(i3,in);
                            end
                        end
                    end
                end
                
            otherwise
                error('Invalid mode')
        end
    else
        error('Invalid dimensions')
    end
end