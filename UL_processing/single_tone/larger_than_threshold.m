function [is_match] = larger_than_threshold(to_be_compare,threshold)
    is_match = 1;
    if to_be_compare < threshold
        is_match = 0;
    end
end