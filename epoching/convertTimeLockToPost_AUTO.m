function [newtrl, react_time] = convertTimeLockToPost_AUTO(ftdata_tmp, trl, tranSamp, plotop, accBadFlag, searchwin, thresh, reactRange, baseMed)
% CONVERTTIMELOCKTOPOST_AUTO detects movement onset from time-locked data automatically.
% This function processes data, applies baseline subtraction, normalizes vectors, and detects movement onset events.
% It also includes options for visualizing the detected movement onset.

%% Initialize Parameters
fsamp = ftdata_tmp.fsample; % Sampling frequency
coder = ftdata_tmp.trial{1}(strncmp(ftdata_tmp.label, 'Coder', 4), :); % Extract coder channel
motInd = find(strncmp(ftdata_tmp.hdr.chantype, 'MotivePos', 4)); % Find motive position indices

% Select motive position data
cfg = [];
cfg.channel = ftdata_tmp.label(motInd);
ftdata_tmp = ft_selectdata(cfg, ftdata_tmp);
XD = ftdata_tmp.trial{1};

% Apply low-pass filter to remove tremor (cutoff frequency = 2 Hz)
fc = 2;
[b, a] = butter(4, fc / (fsamp / 2));
XD = filtfilt(b, a, XD')';

% Calculate distance from baseline and normalize
XD_dist = vecnorm(XD - baseMed);
XD_dist = normalize(XD_dist);

%% Initialize Movement Time and Reaction Time
movtime = [];
react_time = [];

%% Loop Through Each Trial to Determine Movement Onset
for tri = 1:size(trl, 1)
    if (trl(tri, 1) > 0) && (trl(tri, 2) <= size(ftdata_tmp.time{1}, 2))
        % Extract trial data
        mdataA = XD(:, trl(tri, 1):trl(tri, 2));
        mdata = XD_dist(:, trl(tri, 1):trl(tri, 2));
        tvec = linspace(searchwin(1) / fsamp, searchwin(2) / fsamp, numel(mdata));

        % Baseline correction
        mdata = mdata - median(mdata(1:(0.25 * fsamp)));

        % Detect movement based on threshold
        Xinds = SplitVec(find(abs(mdata) > thresh), 'consecutive');
        Xinds = bridgeCellSplits(Xinds, 0.125 * fsamp);
        blipN = 0.125 * fsamp;
        Xinds(cellfun(@length, Xinds) < blipN) = [];

        if ~isempty(Xinds)
            Xinds = Xinds{1};
        end

        offset = 0;

        % Check if movement onset is valid
        if ~isempty(Xinds) && ((Xinds(1) - offset) > 0)
            movlocal = Xinds(1) - offset;

            % Plotting movement detection if visualization is enabled
            if plotop == 1
                subplot(2, 1, 1)
                plot(tvec, mdata)
                hold on; plot([tvec(movlocal), tvec(movlocal)], [min(mdata(:)), max(mdata(:))], 'k--')

                subplot(2, 1, 2)
                plot(tvec, mdataA);
                hold on; plot([tvec(movlocal), tvec(movlocal)], [min(mdataA(:)), max(mdataA(:))], 'k--')
            end

            % Record movement and reaction times
            movtime(tri) = trl(tri, 1) + movlocal;
            react_time(tri) = 1000 * tvec(movlocal);
        else
            movtime(tri) = nan;
            react_time(tri) = nan;
        end

        % Validate reaction time within acceptable range
        if (react_time(tri) < reactRange(1)) || (react_time(tri) > reactRange(2))
            movtime(tri) = nan;
            react_time(tri) = nan;
        end
    else
        warning('Trial definition is longer than recording, replacing with NaNs')
        movtime(tri) = nan;
        react_time(tri) = nan;
    end
end

%% Remove Trials with NaN Movement Times
trl(isnan(movtime), :) = [];
react_time(isnan(movtime)) = [];
movtime(isnan(movtime)) = [];

%% Create New Trial Structure
newtrl = [movtime + tranSamp(1); movtime + tranSamp(2); repmat(tranSamp(1), 1, size(movtime, 2)); trl(:, 4:end)']';

%% Final Plotting if Visualization is Enabled
if plotop == 1
    ax(1) = subplot(2, 1, 1);
    plot(ftdata_tmp.time{1}, ftdata_tmp.trial{1}', 'LineWidth', 2);

    ax(2) = subplot(2, 1, 2);
    X = ftdata_tmp.trial{1} - baseMed;
    mX = [min(X(:)), max(X(:))];
    a = plot(ftdata_tmp.time{1}, X', 'LineWidth', 2);

    yyaxis right
    plot(ftdata_tmp.time{1}, coder', 'LineWidth', 1);

    yyaxis left
    for tri = 1:size(newtrl, 1)
        try
            hold on
            plot([ftdata_tmp.time{1}(newtrl(tri, 1)), ftdata_tmp.time{1}(newtrl(tri, 2))], [mX(1), mX(1)], 'k--')
            plot([ftdata_tmp.time{1}(newtrl(tri, 1)), ftdata_tmp.time{1}(newtrl(tri, 2))], [mX(2), mX(2)], 'k--')
            plot([ftdata_tmp.time{1}(newtrl(tri, 1)), ftdata_tmp.time{1}(newtrl(tri, 1))], mX, 'k--')
            plot([ftdata_tmp.time{1}(newtrl(tri, 1) - tranSamp(1)), ftdata_tmp.time{1}(newtrl(tri, 1) - tranSamp(1))], mX, 'k-')
            plot([ftdata_tmp.time{1}(newtrl(tri, 2)), ftdata_tmp.time{1}(newtrl(tri, 2))], mX, 'k--')
        catch
            disp('Trial definition probably exceeds length of data!')
        end
    end
    legend(a, ftdata_tmp.label)
    xlabel('Time (s)'); ylabel('Movement (mm)')
    linkaxes(ax, 'x')
end

end
