function stat = computeSourceTStatPaired(sourcesaveA,sourcesaveB,fieldn)
% This function computes the TStat of a source space image (repeats saved
% as cells in sourcesave)
stat = sourcesaveA{1};
XTmpA = [];
for i = 1:numel(sourcesaveA)
     XTmpA(:,i) = sourcesaveA{i}.(fieldn)(:);
end

XTmpB = [];
for i = 1:numel(sourcesaveB)
     XTmpB(:,i) = sourcesaveB{i}.(fieldn)(:);
end

% Compute t-statistic
[~,~,~,tstat] = ttest2(XTmpA',XTmpB');
XTstat = tstat.tstat;
XTstat = reshape(XTstat,size(sourcesaveA{1}.(fieldn)));
stat.stat = XTstat;
stat.mask = abs(XTstat)>1.96;