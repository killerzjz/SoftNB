function [wav_sync] = t_sync(wave_in,NRep)
    wave_in = wave_in-mean(real(wave_in))-1i*mean(imag(wave_in));
    num_of_points = 128;
    is_found = 0;
    pattern_length = NRep*960*32;
    for i = 1:NRep*960*32*2
        cor_var1(i) = var(wave_in(i:i+num_of_points));
    end
    max_var1 = max(abs(cor_var1));
    wave_in = wave_in*sqrt(0.25/max_var1);
    max_var1 = 0.25;
    var_threshold = 0.8*mean(max_var1);
    i = 1;
    j = 1;
    startsample = [];
    while i < length(wave_in)-pattern_length-num_of_points
        if larger_than_threshold(abs(var(wave_in(i:i+num_of_points))),var_threshold) && larger_than_threshold(abs(var(wave_in(i+pattern_length-num_of_points:i+pattern_length))),var_threshold) && larger_than_threshold(abs(var(wave_in(i+pattern_length/2-num_of_points:i+pattern_length/2))),var_threshold)
            if i+floor(pattern_length*1.05)-1 < length(wave_in)
                wav_sync_t = wave_in(i:i+floor(pattern_length*1.01)-1);
                is_found = 1;
                i = i + pattern_length - num_of_points;
            else
                wav_sync_t = wave_in(i:end);
                is_found = 1;
                i = i + pattern_length - num_of_points;
            end
            for start_index = 1:num_of_points
                corr(start_index) = 0;
                for lag = 0:9
                    corr(start_index) = corr(start_index)+wav_sync_t(start_index+lag)*conj(wav_sync_t(start_index+lag+128))/abs(wav_sync_t(start_index+lag)*wav_sync_t(start_index+lag+128));
                end
            end
            [~,start] = max(corr);
            start = floor(2*start/7);
            if start == 0
                start = 1;
            end
            startsample(j) = i+start-1-pattern_length+num_of_points;
            wav_sync_t = wav_sync_t(start:start+pattern_length-1);
            wav_sync(j,1:pattern_length) = conj(wav_sync_t');
            j = j + 1;
        end
        i = i + 1;
    end
    if is_found == 0
        wav_sync = [];
    end
end