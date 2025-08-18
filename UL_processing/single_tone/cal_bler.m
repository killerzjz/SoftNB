function[bler] = cal_bler(get_bits,real_bits,NRep)
    size_rx = size(get_bits);
    block_len = 448;
    NErrPackets = 0;
    n_packets = size_rx(1);
    for i_packet = 1:n_packets
        NErrBlocks = 0;
        for i_block = 0:NRep-1
            get_block = get_bits(i_packet,i_block*block_len+1:i_block*block_len+448);
            real_block = real_bits(i_block*block_len+1:i_block*block_len+448);
            if cal_ber(get_block,real_block) > 0
                NErrBlocks = NErrBlocks + 1;
            end
        end
        if NErrBlocks == NRep
            NErrPackets = NErrPackets + 1;
            disp(i_packet);
        end
    end
    bler = NErrPackets/n_packets;
end