function after_rm_cp = rm_cp(msg_in)
    BLOCK_SIZE = 960;        
    HEAD_REMOVE_LEN = 10;    
    TAIL_REMOVE_LEN = 9;    
    REMOVE_PATTERN = [138, 275, 412, 549, 686, 823];
    
    len = length(msg_in);
    remove_mask = false(1, len); 
    
    num_blocks = floor(len / BLOCK_SIZE);
    for block_id = 0:num_blocks-1
        start_idx = block_id * BLOCK_SIZE + 1;
        remove_mask(start_idx : start_idx + HEAD_REMOVE_LEN - 1) = true;
        for offset = REMOVE_PATTERN
            pos = start_idx + offset;
            if pos + TAIL_REMOVE_LEN - 1 <= len
                remove_mask(pos : pos + TAIL_REMOVE_LEN - 1) = true;
            end
        end
    end
    last_start = num_blocks * BLOCK_SIZE + 1;
    if last_start + 138 <= len
        if last_start + HEAD_REMOVE_LEN - 1 <= len
            remove_mask(last_start : last_start + HEAD_REMOVE_LEN - 1) = true;
        end
        last_start = last_start + 138; 
    end

    while last_start + 136 <= len 
        end_idx = min(last_start + TAIL_REMOVE_LEN - 1, len);
        remove_mask(last_start : end_idx) = true;
        last_start = last_start + 137;
    end
    
    after_rm_cp = msg_in(~remove_mask); 
end