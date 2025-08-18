function[symbols] = NB_demod(waveform_in)
symbol_len = floor(length(waveform_in)/128);
for i = 1:symbol_len
    samples = waveform_in((i-1)*128+1:i*128);
    fft_res = fft(samples);
    [~,idx] = max(abs(fft_res));
    symbols(i) = fft_res(idx);
end
end