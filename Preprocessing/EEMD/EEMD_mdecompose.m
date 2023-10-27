function [Y,Imfnumber] = EEMD_mdecompose(X,NoiseLevel,NE,fs,numImf)

[numch,sample] = size(X);
% imfnumber = fix(log2(sample));

Y = [];
for nch = 1 : numch
    disp(['Channel ', num2str(nch)]);
    sigma = NoiseLevel;  
    allmode = rcada_eemd(X(nch,:),sigma,NE,numImf);
    [imfnumber,l] = size(allmode);
    R(nch,:) = X(nch,:)-sum(allmode);  % % residue
    Imfnumber = imfnumber + 1;
    allmode(Imfnumber,:) = R(nch,:);
    X_imfnch = allmode;
    Y = cat(1,Y,X_imfnch);
end

end

