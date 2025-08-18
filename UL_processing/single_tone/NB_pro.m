%Top func
%wave_in_all is signal captured by SDR devices (sample is HackRF One TX and USRP N210 RX, NRep = 4) 
%phaseshift is the phase rotation ruled by the protocols
%NRep is the times of retransmission
function[ber,bler] = NB_pro(wave_in_all,phaseshift,NRep,bits_real)
    wave_sync_t = t_sync(wave_in_all,NRep);
    sizeofwave = size(wave_sync_t);
    n_waves = sizeofwave(1);
    switch NRep
        case 1
            fft_ref = 29277;
        case 2
            fft_ref = 58552;
        case 4
            fft_ref = 117050;
        case 8
            fft_ref = 234465;
        case 16
            fft_ref = 468929;
    end
    for i = 1:n_waves
        wave_one = wave_sync_t(i,1:32*NRep*960);
%         wave_one = t_sync(wave_in);
        sig_len = length(wave_one);
        wave_one = conj(wave_one');
        [~,fft_peak] = max(abs(fft(wave_one)));
        f_offset = fft_ref - fft_peak;
        k = 1:sig_len;
        k = k';
        move_sig = exp(2*pi*1i*f_offset/sig_len*k);
        wave_moved = move_sig.*wave_one;
%        wave_moved = wave_one;
        w_rcp = rm_cp(wave_moved);
        w_dero = derotation(w_rcp,phaseshift);
        w_extended = w_dero;
        % Synchronization with traditional Costas

%         carrierSync = comm.CarrierSynchronizer( ...
%         'SamplesPerSymbol',128, ...
%         'Modulation','QPSK');
%         [w_sync,phaseEstimate] = carrierSync(w_extended);
%         if i== 10
%             save w_sync w_sync;
%             save phaseEstimate phaseEstimate;
%         end

        %Synchronization with local Costas
        w_sync = costas_local(w_extended);
        symbols = NB_demod(w_sync);
        ro1 = (-1-1i)/(symbols(4)/abs(symbols(4)));
        ro2 = (1+1i)/(symbols(11)/abs(symbols(11)));
        ro3 = (-1-1i)/(symbols(18)/abs(symbols(18)));
        ro4 = (-1-1i)/(symbols(25)/abs(symbols(25)));
        ro = (ro1+ro2+ro3+ro4)/4;
        symbols = symbols*ro;
        bits(i,1:448*NRep) = qpskstob(symbols);
    end
    ber = cal_ber(bits,bits_real);
    bler = cal_bler(bits,bits_real,NRep);
end