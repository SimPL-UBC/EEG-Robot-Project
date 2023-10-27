function [reconstructed_clean] = EEMD_CCA(input_eeg, num_IMFs, fs)
    addpath('./EEMD/'); addpath('./CCA/')
    input_eeg = double(input_eeg);
    
    [postEEMD, num_IMFs] = EEMD_mdecompose(input_eeg, 0.2, ...
        10, fs, num_IMFs); 
    postEEMD = single(postEEMD);
    
    [decomposed_data, B_mc, W_mc] = myCCA(postEEMD, fs, 1);
    decomposed_data = real(decomposed_data);
    
    nrows = size(decomposed_data, 1);
    clean_data = zero_artifacts(...
        decomposed_data, round(nrows/2):nrows);
    
    clean_data = inv(W_mc{1,1}')*inv(B_mc(:,:,1))*clean_data;
    reconstructed_clean = real(multiEEG_1_recon(clean_data, clean_data, input_eeg, num_IMFs, fs));
end