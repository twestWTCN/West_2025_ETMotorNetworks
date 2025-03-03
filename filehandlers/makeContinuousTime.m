function ftdata = makeContinuousTime(ftdata)
% This function will ensure that the time vectors for each trial are
% consecutive

tL = cellfun(@length,ftdata.trial); % lengths
ftdata.trial(tL<1) = [];
try
    ftdata.trialinfo(tL<1,:) = [];
catch
%     ftdata.sampleinfo(tL<1,:) = [];
end
fsamp = ftdata.fsample;
initT = 0; T = []; X = [];
for tr = 1:numel(ftdata.trial)
    X(tr) = size(ftdata.trial{tr},2);
    T{tr} = linspace(initT+(1/fsamp), initT+(X(tr)/fsamp),X(tr));
    initT = T{tr}(end);
end

ftdata.time = T;