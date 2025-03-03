function [freq, meanBase] = normalizeSpectrogramByBaseline(baseline, freq, baselinetype, meanBase, trlmeanBase)
%% Function normalizeSpectrogramByBaseline
% This function normalizes the spectrogram data by a specified baseline.
% The baseline normalization method can be of different types, including absolute, relative, relchange, normchange, dB, or zscore.
%
% Arguments:
% baseline - A vector specifying the start and end times for baseline normalization
% freq - Frequency data structure containing power spectrum
% baselinetype - String specifying the type of normalization to be applied
% meanBase - (Optional) Precomputed baseline mean values
% trlmeanBase - (Optional) Flag indicating whether to use the mean over trials (default = 1)
%
% Outputs:
% freq - Frequency data with normalized power spectrum
% meanBase - The mean values used for normalization
if nargin < 5
    trlmeanBase = 1; % This uses the mean over the trials
end

% Adapted from ft_freqbaseline
timeVec = freq.time;
data = freq.powspctrm;

% Calculate baseline mean if not provided
if nargin < 4 || isempty(meanBase)
    baselineTimes = (timeVec >= baseline(1) & timeVec <= baseline(2));
    if numel(size(data)) == 3
        meanBase = nanmean(data(:, :, baselineTimes), 3);
        meanVals = repmat(meanBase, [1 1 size(data, 3)]);
    else % if trials are still present
        meanBase = nanmean(nanmean(data(:, :, :, baselineTimes), 1), 4);
        if trlmeanBase == 0
            meanVals = repmat(nanmean(data(:, :, :, baselineTimes), 4), [1 1 1 size(data, 4)]);
        else
            meanVals = repmat(meanBase, [size(data, 1) 1 1 size(data, 4)]);
        end
    end
else
    disp('Using precomputed baseline!')
    if numel(size(data)) == 3
        meanVals = repmat(meanBase, [1 1 size(data, 3)]);
    else % if trials are still present
        meanVals = repmat(meanBase, [size(data, 1) 1 1 size(data, 4)]);
    end
end

% Apply the baseline normalization
if strcmp(baselinetype, 'absolute')
    data = data - meanVals;
elseif strcmp(baselinetype, 'relative')
    data = data ./ meanVals;
elseif strcmp(baselinetype, 'relchange')
    data = (data - meanVals) ./ meanVals;
elseif strcmp(baselinetype, 'normchange') || strcmp(baselinetype, 'vssum')
    data = (data - meanVals) ./ (data + meanVals);
elseif strcmp(baselinetype, 'db')
    data = 10 * log10(data ./ meanVals);
elseif strcmp(baselinetype, 'zscore')
    stdVals = repmat(nanstd(data(:, :, baselineTimes), 1, 3), [1 1 size(data, 3)]);
    data = (data - meanVals) ./ stdVals;
else
    error('Unsupported method for baseline normalization: %s', baselinetype);
end

% Update the frequency structure with the normalized data
freq.powspctrm = data;
end