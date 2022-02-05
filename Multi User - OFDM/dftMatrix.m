function D = dftMatrix(N,arrayType)

    switch lower(arrayType)
%         case 'square'
%             for m= 0:sqrt(N)-1
%                 for n= 0:sqrt(N)-1
%                     D(:,m+1) = ( 1/sqrt(N) )*exp(-1i*omega*[0:N-1]');
%                 end
%             end
%             y = y.'/sqrt(N);

        case 'linear'
            for m= 0:N-1
                omega = 2*pi*m/N;
%             for m= 0:N-1
%                 omega = 2*pi*m/N;
                D(:,m+1) = ( 1/sqrt(N) )*exp(-1i*omega*[0:N-1]');
            end


        otherwise
            error('Array type not defined!')
            
            
    count = 0;
    
%     for 
%         for phi = 0: 2*pi/Nr : 2*pi - 2*pi/Nr
%             count
%             D(:,1) = array_response(arrayType,AoD_Az(cl,ray),AoD_El(cl,ray),Nt);
%         end 
%     end

end