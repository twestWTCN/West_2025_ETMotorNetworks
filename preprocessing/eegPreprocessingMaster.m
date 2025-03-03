function [data_out, badTr] = eegPreprocessingMaster(R, data_in, plotop)
    % EEGPreprocessingMaster: Preprocesses EEG data including filtering, artifact rejection,
    % ICA, and rereferencing.
    % Parameters:
    %   R - Configuration structure
    %   data_in - FieldTrip data structure containing EEG data
    %   plotop - (Optional) Set to 1 to visualize data; default is 0
    % Returns:
    %   data_out - Processed EEG data
    %   badTr - Trials identified as bad

    if nargin < 3
        plotop = 0;
    end

    %% Step 1: Bugfix for Channel Type Names
    data_in.hdr.chantype(strcmp(data_in.hdr.chantype, 'eeg')) = {'EEG'};

    % Load electrode positions
    load('acticap128_elec', 'elec');
    data_in.elec = elec;

    % Get sensor neighbours
    cfgnbr = [];
    cfgnbr.method = 'template';
    cfgnbr.template = 'acticap128_neighbours.mat';
    neighbours = ft_prepare_neighbours(cfgnbr);

    data_out = data_in;
    eegSel = strncmp(data_in.hdr.chantype, 'EEG', 3);

    %% Step 2: Remove Bad Data Using Visual Inspection
    cfg = [];
    cfg.channel = data_in.label(eegSel(1:5:end));  % Inspect a subset of EEG channels
    cfg = ft_databrowser(cfg, data_in);
    cfg.artfctdef.reject = 'nan';
    data_in = ft_rejectartifact(cfg, data_in);
    [data_in, badGuidePre] = replaceNaNRandom(data_in);

    %% Step 3: Set FCz as Implicit Reference and Update Header
    cfg = [];
    cfg.implicitref = 'FCz';
    data_in = ft_preprocessing(cfg, data_in);
    data_in.hdr.chantype{end+1} = 'EEG';
    data_in.hdr.chanunit{end+1} = 'unknown';
    data_in.hdr.label{end+1} = 'FCz';
    eegSel = strncmp(data_in.hdr.chantype, 'EEG', 3);

    %% Step 4: High-Pass and Low-Pass Filtering
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 2;  % High-pass filter at 2 Hz
    cfg.hpfilttype = 'firws';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 98;  % Low-pass filter at 98 Hz
    cfg.lpfilttype = 'firws';
    cfg.demean = 'no';
    cfg.dftfreq = [50, 100, 150];
    cfg.dftreplace = 'neighbour';
    cfg.dftbandwidth = [1, 1, 1];
    cfg.dftneighbourwidth = [1, 1, 1];
    cfg.channel = data_in.label(eegSel);
    cfg.reref = 'yes';  % Rereference to FCz
    cfg.refchannel = 'FCz';
    data_in = ft_preprocessing(cfg, data_in);

    %% Step 5: Remove the Implicit Reference FCz
    cfg = [];
    cfg.channel = setdiff(data_in.label, 'FCz');
    data_in = ft_selectdata(cfg, data_in);

    %% Step 6: Fix Bad Channels and Mark Bad Trials
    cfg = [];
    cfg.keepchannel = 'repair';
    cfg.neighbours = neighbours;
    cfg.keeptrial = 'zero';
    data_in = ft_rejectvisual(cfg, data_in);

    badTr = cellfun(@(x) all(~any(x(:))), data_in.trial, 'UniformOutput', 0);
    badTr = [badTr{:}];
    cfg = [];
    cfg.trials = find(~badTr);
    data_in = ft_selectdata(cfg, data_in);
    data_out = ft_selectdata(cfg, data_out);  % Update output for good trials
    data_in = correctHeaderLabels(data_in);

    %% Step 7: Remove Remaining Bad Data via Visual Inspection
    eegSel = strncmp(data_in.hdr.chantype, 'EEG', 3);
    cfg = [];
    cfg.channel = data_in.label(eegSel(1:5:end));
    cfg = ft_databrowser(cfg, data_in);
    cfg.artfctdef.reject = 'nan';
    data_in = ft_rejectartifact(cfg, data_in);
    [data_in, badGuidePost] = replaceNaNRandom(data_in);

    %% Step 8: Compute ICA
    dat = cat(2, data_in.trial{:});
    dat(isnan(dat)) = 0;
    n_ic = fix(rank(dat) / 2);
    cfg = [];
    cfg.method = 'runica';
    cfg.numcomponent = n_ic;
    comp = ft_componentanalysis(cfg, data_in);

    % Visualize ICA Components
    cfg = [];
    cfg.layout = 'acticap128';
    cfg.viewmode = 'component';
    cfg.zlim = 'maxmin';
    cfg.compscale = 'global';
    cfg.contournum = 6;
    cfg.artifactalpha = 0.8;
    ft_databrowser(cfg, comp);

    % Reject Specified Components
    badcomps = input('Which components to remove? e.g., [1 2 3]?       ');
    cfg = [];
    cfg.component = badcomps;
    data_in = ft_rejectcomponent(cfg, comp, data_in);

    %% Step 9: Repair Any Remaining Bad Channels
    data_in = repairBadChannel(R, data_in);

    %% Step 10: Apply Common Average Reference
    cfg = [];
    cfg.refmethod = 'avg';
    cfg.refchannel = 'all';
    data_in = ft_preprocessing(cfg, data_in);

    %% Step 11: Replace Bad Data with NaNs
    data_in = putbackNans(data_in, badGuidePre);
    data_in = putbackNans(data_in, badGuidePost);

    %% Step 12: Remove Any Remaining Bad Channels
    cfg = [];
    cfg.channel = [data_in.hdr.label; data_out.hdr.label(~strncmp(data_out.hdr.chantype, 'EEG', 3))];
    data_out = ft_selectdata(cfg, data_out);
    data_out = correctHeaderLabels(data_out);

    eegSel = find(strncmp(data_out.hdr.chantype, 'EEG', 3));
    for tr = 1:numel(data_out.trial)
        data_out.trial{tr}(eegSel, :) = data_in.trial{tr};
    end

    % Final Header Fix
    data_out = correctHeaderLabels(data_out);
end

%% Helper Functions
function data = correctHeaderLabels(data)
    cl_org = data.hdr.label;
    cl_new = data.label;
    [c, ia] = setdiff(cl_org, cl_new);
    data.hdr.label(ia) = [];
    data.hdr.nChans = numel(data.label);
    try
        data.hdr.chantype(ia) = [];
        data.hdr.chanunit(ia) = [];
    catch
        disp('hdr information missing');
    end
end

function [data_in, badGuide] = replaceNaNRandom(data_in)
    % Replace NaNs with random values from non-bad indices
    for tr = 1:numel(data_in.trial)
        X = data_in.trial{tr};
        badInds = find(isnan(X(1, :)));
        badGuide{tr} = zeros(1, size(X, 2));
        badGuide{tr}(badInds) = 1;
        badSec = SplitVec(badInds, 'consecutive');
        for bs = 1:numel(badSec)
            X(:, badSec{bs}) = X(:, getRandInd(size(X, 2), size(badSec{bs}, 2), badInds));
        end
        data_in.trial{tr} = X;
    end
end

function X = getRandInd(dataLength, secLength, badInds)
    % Generate random indices, ensuring no overlap with bad indices
    flag = 0;
    while flag == 0
        X = randi(dataLength);
        X = X:X + (secLength - 1);
        if any(intersect(badInds, X)) || X(end) > dataLength
            flag = 0;
        else
            flag = 1;
        end
    end
end

function data_in = putbackNans(data_in, badGuide)
    % Replace bad segments with NaNs
    for tr = 1:numel(data_in.trial)
        X = data_in.trial{tr};
        X(:, find(badGuide{tr})) = nan;
        data_in.trial{tr} = X;
    end
end

function ftdata = repairBadChannel(R, ftdata)
    badch = find(sum(isnan(ftdata.trial{1}), 2));
    if ~isempty(badch)
        cfg = [];
        cfg.badchannel = ftdata.label(badch);
        cfg.neighbours = R.neighbours;
        ftdata = ft_channelrepair(cfg, ftdata);
    end
end
