function [data_out, badTrials] = artefact_rejection(R, data_in, type)
% ARTEFACT_REJECTION performs artefact rejection on the given data.
% Depending on the type of artefact rejection ('visual' or 'auto'), this function either
% visually assists the user in identifying bad trials or automatically detects artefacts.

% Initialize output
data_out = data_in;

% Determine channel selection based on type
if strcmp(type, 'ALL')
    chSel = 1:numel(data_in.label);
else
    chSel = false(size(data_in.hdr.chantype, 1), 1);
    for i = 1:numel(type)
        chSelInd = strcmp(data_in.hdr.chantype, type{i});
        chSel = chSel | chSelInd;
    end    
end

%% Select the relevant channels
cfg = [];
cfg.channel = data_in.label(chSel);
data_in = ft_selectdata(cfg, data_in);

switch R.artRej.meth
    case 'visual'
        %% Visual Rejection
        cfg = [];
        cfg.keeptrial = 'no';
        switch type{1}
            case 'EEG'
                cfg.keepchannel = 'repair';
                cfg.neighbours = R.neighbours;
            otherwise
                cfg.keepchannel = 'no';
        end
        % Perform visual rejection
        [data_in, badTrials] = ft_rejectvisual(cfg, data_in);
        
    case 'auto'
        %% Automatic Rejection
        % Reject trials based on high variance
        trialLevel = computeARRejectFeats(data_in, 'var');
        trialLevel = mean(trialLevel);
        badTrials{1} = trialLevel > (mean(trialLevel(:)) + (2.33 * std(trialLevel(:))));
        
        % Reject trials based on maximum Z-value
        trialLevel = computeARRejectFeats(data_in, 'maxzvalue');
        channelLevel = max(trialLevel, [], 2);
        badChannels{1} = channelLevel > R.artRej.maxZScoreEpsCh;
        badTrials{2} = max(trialLevel) > R.artRej.maxZScoreEps;
        
        % Reject trials based on flatness
        badTrials{3} = max(trialLevel) < R.artRej.maxZScoreEpsFlat;
        
        % Reject trials based on muscle artefact detection
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [52 98];
        cfg.bpfilttype = 'but';
        data_tmp = ft_preprocessing(cfg, data_in);
        trialLevel = computeARRejectFeats(data_tmp, 'maxzvalue');
        channelLevel = max(trialLevel, [], 2);
        badChannels{2} = channelLevel > R.artRej.maxZScoreMusc;
        badTrials{4} = max(trialLevel) > R.artRej.maxZScoreMusc;
        
        % Combine trial rejection results
        badTrials = badTrials{1} + badTrials{2} + badTrials{3} + badTrials{4};
        badChannels = logical(badChannels{1} + badChannels{2});
        fprintf('Removing %.f bad trials automatically! \n', sum(badTrials));
        
        % Select only good trials
        cfg = [];
        cfg.trials = find(~badTrials);
        data_in = ft_selectdata(cfg, data_in);
        
        % Handle bad channels
        if any(badChannels)
            switch type
                case 'EEG'
                    cfg = [];
                    cfg.badchannel = data_in.label(badChannels);
                    cfg.neighbours = R.neighbours;
                    data_in = ft_channelrepair(cfg, data_in);
                otherwise
                    cfg = [];
                    cfg.channel = data_in.label(~badChannels);
                    cfg.keepchannel = 'no';
                    data_in = ft_selectdata(cfg, data_in);
            end
        end
        
        % Ensure badTrials is an index vector
        badTrials = find(badTrials);
        error('Check this autorej for badTr compat (must output index not logical)');
end

%% Recombine Data by Selecting Trials from Original Data
cfg = [];
cfg.channel = union(data_in.label, data_out.label(~chSel)); % Select identified bad channels
cfg.trials = setdiff(1:numel(data_out.trial), badTrials); % Select trials not marked as bad
data_out = ft_selectdata(cfg, data_out);

% Correct the header to account for removed channels
data_out = correctHeaderLabels(data_out);

end

function data = correctHeaderLabels(data)
% CORRECTHEADERLABELS updates the header labels to reflect removed channels.
cl_org = data.hdr.label;
cl_new = data.label;

[c, ia] = setdiff(cl_org, cl_new);
data.hdr.label(ia) = [];
data.hdr.nChans = numel(data.label);
try
    data.hdr.chantype(ia) = [];
    data.hdr.chanunit(ia) = [];
catch
    disp('Header information missing');
end

end