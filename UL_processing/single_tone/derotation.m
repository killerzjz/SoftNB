function[after_derotation] = derotation(after_rmcp,phaseshift)
symbol_len = length(after_rmcp)/128;
for i = 1:symbol_len
    ro_wave = complex(zeros(128,1));
    ro_wave(123) = phaseshift(i);
    ro_wave = ifft(ro_wave);
    after_derotation((i-1)*128+1:i*128) = after_rmcp((i-1)*128+1:i*128)./ro_wave;
end
    after_derotation = conj(after_derotation');
end