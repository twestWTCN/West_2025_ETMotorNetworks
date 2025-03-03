function [peakPow, peakFrq, tremPow, tremVar, tmpFx, accData, trem_Env] = computeTremorProperties(ftdata_cat, part, keeptrials)
% Compute tremor properties from accelerometer data
% Inputs:
%   ftdata_cat - FieldTrip data structure containing the signals
%   part - specific parts of the task where tremor occurs (optional)
%   keeptrials - flag to keep individual trial data (optional, default = 1)
% Outputs:
%   peakPow - peak power within the tremor frequency range
%   peakFrq - peak frequency within the tremor frequency range
%   tremPow - RMS power of the tremor signal
%   tremVar - variance of the tremor envelope
%   tmpFx - frequency analysis results
%   accData - processed accelerometer data
%   trem_Env - envelope of the tremor signal

if nargin < 2
    part = [];
end
if nargin < 3
    keeptrials = 1;
end

% Select the accelerometer data
cfg = [];
cfg.channel = {'AccX', 'AccY', 'AccZ'};
accData = ft_selectdata(cfg, ftdata_cat);

% Select the parts of the task where tremor occurs
if numel(part) == 1
    cfg = [];
    cfg.latency = getTremorParts(part);
    accData = ft_selectdata(cfg, accData);
elif numel(part) == 3
    cfg = [];
    cfg.latency = [part(1), part(2)];
    cfg.trials = part(3);
    accData = ft_selectdata(cfg, accData);
end

% Convert signals from volts to acceleration and detrend
accData.trial = cellfun(@(x) signal2accel(x), accData.trial, 'UniformOutput', false);
accData.trial = cellfun(@(x) detrend(x', 3)', accData.trial, 'UniformOutput', false);

% Perform PCA on the acceleration data
[~, accData] = getFTAccPca(accData, [2 12]);

% Select only the first principal component
cfg = [];
cfg.channel = 'pca001';
accData = ft_selectdata(cfg, accData);

% Compute the power spectrum
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.foi = 2:0.5:24;
cfg.tapsmofrq = 1.5;
if keeptrials
    cfg.keeptrials = 'yes';
else
    cfg.keeptrials = 'no';
end
cfg.pad = 'nextpow2';
cfg.channel = 'pca001';
tmpFx = ft_freqanalysis(cfg, accData);

% Extract peak frequency and power from the power spectrum
[peakFrq, peakPow] = powTremorProc(tmpFx, keeptrials);

% Compute tremor properties such as RMS power, envelope, and variance
[~, ~, tremPow, trem_Env, tremVar] = bandpass_and_rms(accData.trial, accData.fsample, mean(peakFrq), 3);

end
