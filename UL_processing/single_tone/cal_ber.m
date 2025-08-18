function[ber] = cal_ber(get_bits,real_bits)
    n_packets = size(get_bits);
    n_packets = n_packets(1);
    for i = 1:n_packets
        err = get_bits(i,1:length(real_bits)) - real_bits;
        ber(i) = sum(err.*err)/length(err);
    end
end