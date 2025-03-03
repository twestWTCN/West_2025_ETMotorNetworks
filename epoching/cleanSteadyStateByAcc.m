function data_out = cleanSteadyStateByAcc(data_in, epsprc)
% cleanSteadyStateByAcc - Filters accelerometer data to extract steady-state segments.
% This function processes accelerometer data to identify and retain segments where
% the signal is relatively steady, using a band-stop filter and statistical measures.
% It removes noisy segments to ensure more stable data for further analysis.

% Identify indices of channels labeled as Accelerometer
accInd = find(strcmp(data_in.hdr.chantype, 'Accelerometer'));

% Apply band-stop filter to remove oscillations in the frequency range [2-15] Hz
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [2 15];
cfg.channel = data_in.label(accInd);
ftdata_block_filt = ft_preprocessing(cfg, data_in);

% Iterate over each trial to identify steady-state segments
for tr = 1:numel(ftdata_block_filt.trial)
    T = ftdata_block_filt.time{tr}; % Time vector for the trial
    X = sum(ftdata_block_filt.trial{tr}); % Sum of accelerometer channels
    X = (X - mean(X)) ./ std(X); % Z-score normalization
    X = smooth(X, 0.05 * ftdata_block_filt.fsample)'; % Smooth data with a window of 5% of sampling rate
    
    % Split the data into segments for further analysis
    N = fix(0.25 * ftdata_block_filt.fsample); % Segment length of 0.25 seconds
    if numel(X) > 2 * N
        K = X(1:end - rem(numel(X), N)); % Truncate to make length divisible by N
        XN = reshape(K, N, []); % Reshape into segments of length N
        
        % Calculate threshold for steady-state identification
        eps = prctile(max(abs(XN)), epsprc); % Upper percentile threshold based on epsprc
        
        % Identify and bridge small gaps between steady-state segments
        Xinds = SplitVec(find(abs(X) < eps), 'consecutive');
        Xinds = bridgeCellSplits(Xinds, 0.5 * ftdata_block_filt.fsample);
                
        % Find the longest steady-state segment
        [~, a] = max(cellfun(@length, Xinds));
        steadyBlock = Xinds{a}(1):Xinds{a}(end); % Continuous steady-state block
        
        % Mark the start and end of the steady-state segment
        plot([T(steadyBlock(1)) T(steadyBlock(1))], [-5 5], 'r--');
        plot([T(steadyBlock(end)) T(steadyBlock(end))], [-5 5], 'r--');
        
        % Save the steady-state segment definition
        trldef(tr, :) = [steadyBlock(1), steadyBlock(end), length(steadyBlock)];
    else
        % If data is too short, use the entire segment
        trldef(tr, :) = [1, numel(X), length(X)];
    end
    clf; % Clear figure for the next trial
end

% Redefine trials based on the identified steady-state segments
cfg = [];
cfg.trl = trldef;
data_out = ft_redefinetrial(cfg, data_in);

% Plot the resulting steady-state segments for verification
for tr = 1:numel(data_out.trial)
    T = data_out.time{tr};
    X = sum(data_out.trial{tr});
    plot(T, X);
    hold on;
end

% Preserve the history field from the input data
data_out.hdr.history = data_in.hdr.history;
clf; % Clear the figure for final cleanup
end
