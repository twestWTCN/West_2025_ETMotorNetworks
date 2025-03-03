function ftdata = checkConsistentDataLength(ftdata,stdthresh)
SL = std(cellfun(@(x) size(x,2),ftdata.trial));
if SL>stdthresh
    warning('Trials in this data have a mismatch in data length! Press any key to continue')
    % yn = input('Do you want to correct the lengths? y/n','s');
    % if yn == 'y'
        len = cellfun(@length,ftdata.trial);
        if min(len)<512
            warning('Length is really short!')
            pause
        end
        ftdata = makeSameLength(ftdata,min(len));
        checkConsistentDataLength(ftdata,0.2)
    % end
end