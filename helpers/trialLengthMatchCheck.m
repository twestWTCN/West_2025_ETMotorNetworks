function trialLengthMatchCheck(FTdata1,FTdata2)

L1 = cellfun(@(x) size(x,2),FTdata1.trial);
L2 = cellfun(@(x) size(x,2),FTdata2.trial);

if numel(L1) ~= numel(L2)
    warning('Data have different trial numbers!')
end

[cellLem,cellInd] = setdiff(L1,L2);

if any(cellLem)
    warning('Your data lengths are not matching!')
end