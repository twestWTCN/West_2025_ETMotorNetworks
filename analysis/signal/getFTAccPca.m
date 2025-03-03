function [ftdata, AccDat] = getFTAccPca(ftdata, tremfreq, filtall)
% Apply PCA to accelerometer data in FieldTrip data structure
% Inputs:
%   ftdata - FieldTrip data structure containing the signals
%   tremfreq - (optional) bandpass filter frequencies [low high] for preprocessing
%   filtall - (optional) flag to filter all components
% Outputs:
%   ftdata - FieldTrip data structure appended with PCA components of accelerometer data
%   AccDat - Preprocessed accelerometer data after PCA

% Select accelerometer channels
cfg = [];
cfg.channel = {'AccX', 'AccY', 'AccZ'};
AccDat = ft_selectdata(cfg, ftdata);

% If tremor frequency range is provided, apply bandpass filter
if nargin > 1
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = tremfreq;
    AccDat = ft_preprocessing(cfg, AccDat);
end

% Remove unnecessary fields if present
if isfield(AccDat, 'grad')
    AccDat = rmfield(AccDat, 'grad');
end
if isfield(AccDat, 'elec')
    AccDat = rmfield(AccDat, 'elec');
end

% Apply PCA to the accelerometer data
cfg = [];
cfg.method = 'pca';
AccDat = ft_componentanalysis(cfg, AccDat);

% Remove components not needed
AccDat = rmfield(AccDat, 'topo');
AccDat = rmfield(AccDat, 'topolabel');
AccDat = rmfield(AccDat, 'unmixing');

% Append PCA components back to the original FieldTrip data
cfg = [];
ftdata = ft_appenddata(cfg, ftdata, AccDat);

end