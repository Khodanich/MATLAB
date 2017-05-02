function z = Correlation(x, y) %% x and y are massives of complex numbers
if( length(x) < length(y) )
    len = length(y); 
    for inc = length(x)+1:length(y)
        x(inc) = 0;
    end
end

if( length(y) < length(x) )
    len = length(x);
    y(length(y)+1:length(x)) = 0;
end

if ( length(x) == length(y) ) 
    len = length(x);
end
    
z_fft = zeros(len, 1);
z     = zeros(len, 1);
    x_fft = fft(x);
    y_fft = fft(y);
   
    z_fft(:, 1) = x_fft(:).*conj(y_fft(:));
   
    zet = ifft(z_fft);
%     z(:, 1) = abs(zet(:, 1));
    z(:, 1) = zet(:, 1);
end
    