function stat = computeSourceTStat(sourcesave,fieldn)
% This function computes the TStat of a source space image (repeats saved
% as cells in sourcesave)
stat = sourcesave{1};
XTmp = [];
for i = 1:numel(sourcesave)
    XTmp(:,i) = sourcesave{i}.(fieldn)(:);
end
% OR Compute t-statistic
XBar = mean(XTmp');
XSE = std(XTmp')./sqrt(size(XTmp,2));
XTstat = XBar./XSE;
XTstat = reshape(XTstat,size(sourcesave{1}.(fieldn)));
stat.stat = XTstat;
stat.mask = abs(XTstat)>1.96;