function [ber,bler]  = test_with_NRep4()
    load("bits_real\bits4.mat");
    load("phaseshift\\phaseshift_NRep4.mat");
    load("sample_data\NRep4.mat")
    [ber,bler] = NB_pro(NRep4,phaseshift_NRep4,4,bits4);
end
