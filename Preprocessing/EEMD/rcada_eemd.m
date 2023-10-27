function [allmode] = rcada_eemd(Y, NoiseLevel, NE, numImf, varargin)
allmode = [];
[Y,NoiseLevel,NE,numImf,runCEEMD,maxSift,typeSpline,toModifyBC,randType,seedNo,IsInputOkay] = parse_checkProperty(Y, NoiseLevel, NE, numImf, varargin);

if(~IsInputOkay)
    fprintf('ERROR : The process is not executed.\n');
    return;
end

if (NoiseLevel == 0)
    allmode = rcada_emd(Y, toModifyBC, typeSpline, numImf, maxSift);
    allmode = allmode';
    return;
end

xsize = size(Y,2);
Ystd = std(Y);	

allmode = zeros(xsize,numImf);

savedState = set_seed(seedNo);	
if (runCEEMD)
   NE = 2*NE; % YHW0202_2011: flip noise to balance the perturbed noise
end

for iii=1:NE  % ensemble loop    
    if (runCEEMD)
        if (mod(iii,2) ~= 0)
            if (randType == 1) % White Noise
               temp = ((2*rand(1,xsize)-1)*NoiseLevel).*Ystd; 
            elseif (randType == 2) % Gaussian Noise
               temp = (randn(1,xsize)*NoiseLevel).*Ystd; 
            end
         else % Even number
           temp = -temp;
        end
    else % runCEEMD = 0
        if (randType == 1)
            temp = (2*rand(1,xsize)-1)*NoiseLevel.*Ystd; % temp is Ystd*[0 1]	 
        elseif (randType == 2)
            temp = randn(1,xsize)*NoiseLevel.*Ystd; % temp is Ystd*[0 1]
        end 
    end
    xend =  Y + temp;
    imf = rcada_emd(xend, toModifyBC, typeSpline, numImf, maxSift);
    allmode = allmode + imf; 
end % iii: ensemble loop

return_seed(savedState);
allmode = allmode/NE;

allmode = allmode'; % 0318_2014
% return; % end eemd
end

function savedState = set_seed(seedNo)
% defaultStream = RandStream.getDefaultStream;
defaultStream = RandStream.getGlobalStream;
savedState = defaultStream.State;
rand('seed',seedNo);
randn('seed',seedNo);

end

function return_seed(savedState)
RandStream.getDefaultStream.State = savedState;
end

function [Y, NoiseLevel, NE, numImf, runCEEMD, maxSift, typeSpline,toModifyBC,randType,seedNo, IsInputOkay] = parse_checkProperty(Y, NoiseLevel, NE, numImf, varargin)
% Default Parameters
runCEEMD = 0; % Original EEMD
maxSift = 10; % maxSift = 10
typeSpline = 2;
toModifyBC = 1;
randType = 2;
seedNo = 1; % now
checkSignal = 0;
IsInputOkay = true;

if(~isempty(varargin{1}))

for iArg = 1 : length(varargin{1});
    
if(iArg == 1)
   runCEEMD = varargin{1}{iArg};
   if(runCEEMD ~= 0 && runCEEMD ~= 1)
    fprintf('ERROR : runCEEMD must be 0 (Off) or 1 (On).\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 2)
   maxSift = varargin{1}{iArg};
   if(maxSift < 1 || (mod(maxSift, 1) ~= 0))
    fprintf('ERROR : Number of Iteration must be an integer more than 0.\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 3)
    typeSpline = varargin{1}{iArg};
    if(typeSpline ~= 1 && typeSpline ~= 2 && typeSpline ~= 3)
    fprintf('ERROR : typeSpline must be 1 (clamped spline); 2 (not a knot spline).\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 4)
    toModifyBC = varargin{1}{iArg};
    if(toModifyBC ~= 0 && toModifyBC ~= 1 && toModifyBC ~= 2)
    fprintf('ERROR : toModifyBC must be 0 (None) ; 1 (modified linear extrapolation); 2 (Mirror Boundary)\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 5)
    randType = varargin{1}{iArg};
    if(randType ~= 1 && randType ~= 2)
    fprintf('ERROR : randType must be 1 (uniformly distributed white noise) ; 2 (gaussian white noise).\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 6)
    seedNo = varargin{1}{iArg};
    if(seedNo < 0 || seedNo >  2^32-1 || (mod(seedNo, 1) ~= 0))
    fprintf('ERROR : The value of seed must be an integer between 0 and 2^32 - 1. \n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 7)
    checkSignal = varargin{1}{iArg};
    if(checkSignal ~= 0 && checkSignal ~= 1)
    fprintf('ERROR : Number of checksignal must be 1 (Yes) or 0 (No).\n');
    IsInputOkay = false;
    return;
    end
end

end

end

if(NoiseLevel == 0)
    fprintf('If NoiseLevel is ZERO, EEMD algorithm will be changed to EMD algorithm.\n');
end
if ((NE < 1) || (mod(NE, 1) ~= 0))
    fprintf('ERROR : Number of Ensemble must be integer more than 0.\n');
    IsInputOkay = false;
    return;
end


[m,n] = size(Y);
if(m ~= 1)
    if((n ~= 1))
       fprintf('ERROR : EMD could not input matrix array !\n');
       IsInputOkay = false;
       return;
    else
        Y =Y';
        xsize = m;
    end
else
    xsize = n;
end

if (checkSignal == 1)
    if((any(isinf(Y(:)) == 1)) || (any(isnan(Y(:)) == 1)))
        fprintf('ERROR : The input signal has NaN or Infinity elements.\n');
        IsInputOkay = false;
        return;
    end
end

if(mod(numImf, 1) ~= 0)
    fprintf('ERROR : numImf must be an integer more than 0. \n');
    IsInputOkay = false;
    return;
end

if (numImf <= 0) % automatic estimating number of imf 
    numImf=fix(log2(xsize));
end

end
