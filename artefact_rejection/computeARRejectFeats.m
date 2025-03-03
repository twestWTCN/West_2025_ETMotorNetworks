function level = computeARRejectFeats(ftdata,metric)
% this script is adapted from Fieldtrip's "rejectvisual_summary.m" script
nchan   = numel(ftdata.label);
ntrl   = numel(ftdata.trial);

if strcmp(metric, 'zvalue') || strcmp(metric, 'maxzvalue')
    % cellmean and cellstd (see ft_denoise_pca) would work instead of for-loops, but they are too memory-intensive
    runsum = zeros(nchan, 1);
    runss  = zeros(nchan, 1);
    runnum = 0;
    for i=1:ntrl
        dat = ftdata.trial{i};
        runsum = runsum + nansum(dat, 2);
        runss  = runss  + nansum(dat.^2, 2);
        runnum = runnum + sum(isfinite(dat), 2);
    end
    mval = runsum./runnum;
    sd   = sqrt(runss./runnum - (runsum./runnum).^2);
end
level = [];
for i=1:ntrl
%     ft_progress(i/ntrl, 'computing metric %d of %d\n', i, ntrl);
    dat = ftdata.trial{i};
    switch metric
        case 'var'
            level(:, i) = nanstd(dat, [], 2).^2;
        case 'min'
            level(:, i) = nanmin(dat, [], 2);
        case 'max'
            level(:, i) = nanmax(dat, [], 2);
        case 'maxabs'
            level(:, i) = nanmax(abs(dat), [], 2);
        case 'range'
            level(:, i) = nanmax(dat, [], 2) - nanmin(dat, [], 2);
        case 'kurtosis'
            level(:, i) = kurtosis(dat, [], 2);
        case '1/var'
            level(:, i) = 1./(nanstd(dat, [], 2).^2);
        case 'zvalue'
            level(:, i) = nanmean( (dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , 2);
        case 'maxzvalue'
            level(:, i) = nanmax( ( dat-repmat(mval, 1, size(dat, 2)) )./repmat(sd, 1, size(dat, 2)) , [], 2);
        otherwise
            ft_error('unsupported method');
    end
end