function [Y] = multiEEG_1_recon(X,X_imfchoosen,X_original,componumber,fs)
% program  by  Xueyuan Xu and Luchang Li.
% multiEEG_recon generates a denoised EEG matrix Y from the X.
% Note: > X_orignal: the original observation matrix of dimension [m by T].
%       > m: sensor number.
%       > T: signal length.
%       > fs: sampling frequency.
%       > X_imfchoosen: the IMF matrix after previous choosen fuction
%       > X: the Component matrix after previous denoised fuction
%       > componumber: number of the component processed of every channel of the original signal X. 
    %       > Y: the denoised EEG signals.
[numch sample] = size(X_original);

for modulenumber = 1: numch
     firstcomponent = 1 + componumber*(modulenumber-1);
     endcomponent = componumber*modulenumber;
     Y(modulenumber,:) = sum(X_imfchoosen(firstcomponent:endcomponent,:));
end 
end

