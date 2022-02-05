function y = array_response(arrayType,a1,a2,N)
    
    switch lower(arrayType)
        case 'square'
            for m= 0:sqrt(N)-1
                for n= 0:sqrt(N)-1
                    y(m*(sqrt(N))+n+1) = exp(1i*pi*(m*sin(a1)*sin(a2) + n*cos(a2)));
                end
            end
            y = y.'/sqrt(N);

        case 'linear'
            for m= 0:N-1
                y(m+1) = exp(1i*pi*(m*sin(a1)));
            end
            y = y.'/sqrt(N);

        otherwise
            error('Array type not defined!')
                
                
%     if arrayType == 'square'
%         for m= 0:sqrt(N)-1
%             for n= 0:sqrt(N)-1
%                 y(m*(sqrt(N))+n+1) = exp( 1i* pi* ( m*sin(a1)*sin(a2) + n*cos(a2) ) );
%             end
%         end
%         y = y.'/sqrt(N);
%     else
%         error('Array type not defined!');
%     end
    
    
%         At(k,:,:) = (1/sqrt(Nt))*exp(sqrt(-1)*[0:Nt-1]'*pi*sin(vec(AoD)'));
%         Ar(k,:,:) = (1/sqrt(Nr))*exp(sqrt(-1)*[0:Nr-1]'*pi*sin(vec(AoA)'));

end