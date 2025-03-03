function [X_filtered, X_rms, X_mean, X_Env, X_std] = bandpass_and_rms(X, fs, peak_freq, bandwidth)
% Bandpass filter the data and compute RMS and related metrics
% Inputs:
%   X - cell array of trials (each cell contains a matrix of signals)
%   fs - sampling frequency in Hz
%   peak_freq - center frequency of interest
%   bandwidth - bandwidth around the center frequency for filtering
% Outputs:
%   X_filtered - filtered data (cell array)
%   X_rms - RMS value of each trial
%   X_mean - mean value of the envelope of each trial
%   X_Env - envelope of each trial
%   X_std - standard deviation of the envelope of each trial

% Detrend the data for each trial
X = cellfun(@(x) detrend(x, 3), X, 'UniformOutput', false);

% Concatenate all trials along the time dimension
concatenated_data = cell2mat(X);
concatenated_data = demean(concatenated_data, 2);

% Define the bandpass FIR filter
low_cutoff = peak_freq - bandwidth;
if low_cutoff < 1
    low_cutoff = 1;
    warning('Low cutoff is <1 Hz, truncating to 1 Hz');
end
high_cutoff = peak_freq + bandwidth;
filter_order = 300; % Adjust based on specific requirements
b = fir1(filter_order, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');

% Apply the bandpass filter using filtfilt for zero-phase filtering
concatenated_filtered = filtfilt(b, 1, concatenated_data')';

% Split the concatenated filtered data back into trials
trial_lengths = cellfun(@(x) size(x, 2), X);
X_filtered = mat2cell(concatenated_filtered, size(concatenated_filtered, 1), trial_lengths);

% Initialize RMS, mean, and std matrices
X_rms = nan(size(X{1}, 1), length(X));
X_mean = nan(size(X{1}, 1), length(X));
X_std = nan(size(X{1}, 1), length(X));

% Compute the RMS, mean, std, and envelope for each trial
X_envelope = [];
for i = 1:length(X)
    % Compute the amplitude envelope
    X_envelope(:, i) = abs(hilbert(X_filtered{i}')');

    % Compute RMS, mean, and standard deviation
    X_rms(:, i) = sqrt(mean(X_filtered{i}.^2, 2));
    X_mean(:, i) = mean(X_envelope(:, i));
    X_std(:, i) = std(X_envelope(:, i));
end

% Convert the envelope data back to cell format
X_Env = mat2cell(X_envelope(:)', size(X_envelope(:)', 1), trial_lengths);

end
