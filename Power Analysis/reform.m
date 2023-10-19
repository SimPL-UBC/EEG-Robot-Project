function output = reform(input)
    ts = linspace(1, length(input), length(input));
    tsOldPre = ts(ts<=4);
    tsOldMid = ts(logical((ts>4).*(ts<=24)));
    tsOldPos = ts(ts>24);
    tsNew = linspace(1, 28, 28);
    tsNewPre = tsNew(tsNew<=4);
    tsNewMid = tsNew(logical((tsNew>4).*(tsNew<=24)));
    tsNewPos = tsNew(tsNew>24);
    pre = input(ts<=4);
    mid = input(logical((ts>4).*(ts<=24)));
    pos = input(ts>24); 
    tempPre = interp1(tsOldPre, pre, tsNewPre);
    tempMid = interp1(tsOldMid, mid, tsNewMid);
    tempPos = interp1(tsOldPos, pos, tsNewPos);
    temp = [tempPre tempMid tempPos];
    output = fillmissing(temp, 'nearest');
end