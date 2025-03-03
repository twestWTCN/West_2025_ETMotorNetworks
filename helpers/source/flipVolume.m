function statFlip = flipVolume(LR,stat,dim)
if strcmp(LR,'R') || strcmp(LR,'?')
    % do nothing
    statFlip = stat;
elseif strcmp(LR,'L')
    % Flip L/R
    XStat = reshape(stat,dim);
    XStatFlip = flip(XStat,1);
    statFlip = XStatFlip(:);
end