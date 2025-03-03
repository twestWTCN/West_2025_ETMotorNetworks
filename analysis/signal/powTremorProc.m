function [peakFrq, peakPow, fxOrg, hz_org, axSel] = powTremorProc(freq, keeptrials)
% Process power spectra to determine tremor properties
% Inputs:
%   freq - frequency domain data (FieldTrip structure)
%   keeptrials - flag indicating whether to keep individual trial data
% Outputs:
%   peakFrq - peak frequency within the tremor band (4-12 Hz)
%   peakPow - peak power within the tremor band (4-12 Hz)
%   fxOrg - original power spectrum data
%   hz_org - frequency vector
%   axSel - channel with maximum power

% Initialize outputs
peakPow = [];
peakFrq = [];

% Loop through each channel to calculate the power spectrum
for i = 1:numel(freq.label)
    hz_org = freq.freq;
    switch freq.dimord
        case 'rpt_chan_freq'
            fx_org{i} = squeeze(nanmean(freq.powspctrm(:, i, :), 1));
        case 'chan_freq'
            fx_org{i} = freq.powspctrm(i, :);
        otherwise
            error('Frequency data does not have the expected dimord!');
    end
    
    % Find peak power and frequency within the tremor range (4-12 Hz)
    peakPow(i) = max(fx_org{i}(hz_org >= 4 & hz_org <= 12));
    [~, loc] = max(fx_org{i}(hz_org >= 4 & hz_org <= 12));
    hz_red = hz_org(hz_org >= 4 & hz_org <= 12);
    peakFrq(i) = hz_red(loc);
end

% Select the channel with maximum peak power
[~, axSel] = max(peakPow);

if keeptrials
    % Keep individual trial data for the selected channel
    fx_org = squeeze(freq.powspctrm(:, axSel, :));
    peakPow = [];
    peakFrq = [];
    
    % Loop through each trial to calculate peak power and frequency
    for tr = 1:size(fx_org, 1)
        fx_org_tr = fx_org(tr, :);
        if any(fx_org_tr)
            peakPow(tr) = max(fx_org_tr(hz_org >= 4 & hz_org <= 12));
            hz_red = hz_org(hz_org >= 4 & hz_org <= 12);
            [~, loc] = max(fx_org_tr(hz_org >= 4 & hz_org <= 12));
            peakFrq(tr) = hz_red(loc);
        else
            peakFrq(tr) = nan;
            peakPow(tr) = nan;
        end
    end
    fxOrg = fx_org;
else
    % Use collapsed data for the selected channel
    peakFrq = peakFrq(axSel);
    peakPow = peakPow(axSel);
    fxOrg = fx_org{axSel};
end

end
