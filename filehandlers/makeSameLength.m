function ftdata = makeSameLength(ftdata,len)
for tr = 1:numel(ftdata.trial)
    ftdata.trial{tr} = ftdata.trial{tr}(:,1:len);
    ftdata.time{tr} = ftdata.time{tr}(1:len);
end
ftdata = rmfield(ftdata,'sampleinfo');
