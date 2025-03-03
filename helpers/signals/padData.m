function ftdata = padData(ftdata,dSize,timer,padval)
X = padarray(ftdata.trial{1},[0 dSize],inf,'pre');
if isnumeric(padval)
    X(X==inf) = padval;
elseif strcmp(padval,'randn')
    X(X==inf) = randn;
end
ftdata.trial{1} = X;
% Construct new time vector
ftdata.time = timer;
ftdata.sampleinfo = [0 numel(ftdata.time{1})];