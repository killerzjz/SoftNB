function[bits] = qpskstob(symbols)
s_len = length(symbols);
for i = 1:s_len
    if real(symbols(i))>0 && imag(symbols(i))>0
        bits(2*i-1) = 1;
        bits(2*i) = 1;
    end
    if real(symbols(i))>0 && imag(symbols(i))<0
        bits(2*i-1) = 1;
        bits(2*i) = 0;
    end
    if real(symbols(i))<0 && imag(symbols(i))>0
        bits(2*i-1) = 0;
        bits(2*i) = 1;
    end
    if real(symbols(i))<0 && imag(symbols(i))<0
        bits(2*i-1) = 0;
        bits(2*i) = 0;
    end
bits = bits';
end