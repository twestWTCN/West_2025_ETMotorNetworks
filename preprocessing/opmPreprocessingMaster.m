function [data_out, badTr] = opmPreprocessingMaster(R, data_in, plotop)
    % OPMPreprocessingMaster: Preprocesses OPM (Optically Pumped Magnetometer) data.
    % Applies filtering, ICA, and AMM, among other preprocessing steps.
    % Parameters:
    %   R - Configuration structure
    %   data_in - FieldTrip data structure containing OPM data
    %   plotop - (Optional) Set to 1 to visualize data; default is 0
    % Returns:
    %   data_out - Processed OPM data
    %   badTr - Trials identified as bad

    if nargin < 3
        plotop = 0;
    end
    
    % Initialize output as input
    data_out = data_in;
    opmSel = strncmp(data_in.label, 'G', 1);
    close all;

    %% Step 1: Inspect Initial Spectra (Optional Plotting)
    if plotop
        cfg = [];
        cfg.channel = data_in.label(opmSel);
        cfg.trial_length = 3;
        cfg.method = 'tim';
        cfg.foi = [0.5, 98];
        cfg.plot = 'yes';
        ft_opm_psd(cfg, data_in);
        title('Pre HFC');
        axis square;
        set(gcf, 'Position', [100.0000, -9.4000, 799.4000, 792.0000]);
        ylim([1e-1, 1e5]);
    end

    %% Step 2: Remove Bad Channels with Bandpass Filtering
    cfg = [];
    cfg.preproc.bpfilter = 'yes';
    cfg.preproc.bpfilttype = 'but';
    cfg.preproc.bpfreq = [1, 98];
    cfg.keepchannel = 'no';
    cfg.keeptrial = 'zero';
    cfg.channel = data_in.label(opmSel);
    data_in = ft_rejectvisual(cfg, data_in);
    trialLengthMatchCheck(data_in, data_out);

    badTr = cellfun(@(x) all(~any(x(:))), data_in.trial, 'UniformOutput', 0);
    badTr = [badTr{:}];

    %% Step 3: Select Good Trials and Correct Channel Labels
    cfg = [];
    cfg.trials = find(~badTr);
    data_in = ft_selectdata(cfg, data_in);
    cfg.channel = union(data_out.label(~opmSel), data_in.label(opmSel));
    data_out = ft_selectdata(cfg, data_out);
    data_in = correctHeaderLabels(data_in);
    trialLengthMatchCheck(data_in, data_out);

    %% Step 4: Remove Remaining Bad Data via Artifact Rejection
    cfg = [];
    cfg.preproc.hpfilter = 'yes';
    cfg.preproc.hpfilttype = 'but';
    cfg.preproc.hpfreq = 2;
    cfg.channel = data_in.label(opmSel(1:2:end));
    cfg = ft_databrowser(cfg, data_in);
    cfg.artfctdef.reject = 'nan';
    data_in = ft_rejectartifact(cfg, data_in);
    [data_in, badGuidePost] = replaceNaNRandom(data_in);
    trialLengthMatchCheck(data_in, data_out);

    %% Step 5: High-Pass Filtering at 1 Hz, Low-Pass Filtering at 98 Hz
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 1;
    cfg.hpfilttype = 'firws';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 98;
    cfg.lpfilttype = 'firws';
    cfg.demean = 'yes';
    cfg.dftfreq = [50, 100, 150];
    cfg.dftreplace = 'neighbour';
    cfg.dftbandwidth = [1, 1, 1];
    cfg.dftneighbourwidth = [1, 1, 1];
    cfg.channel = data_in.label(opmSel);
    data_in = ft_preprocessing(cfg, data_in);

    %% Step 6: Fix Gradiometer Definitions
    data_in.grad.chantype = repmat({data_in.grad.chantype}, 1, size(data_in.grad.chanpos, 1))';
    data_in.grad.chanunit = repmat({'T'}, 1, size(data_in.grad.chanpos, 1))';

    %% Step 7: Run ICA and Remove Artifacts
    comp = runICA(data_in);
    badcomps = input('Which ICA components to remove? e.g., [1 2 3]?       ');
    trialLengthMatchCheck(data_in, data_out);

    cfg = [];
    cfg.component = badcomps;
    data_in = ft_rejectcomponent(cfg, comp, data_in);

    %% Step 8:AMM denoising Using SPM
    tmpData = data_in;
    tmpData.trial = {[data_in.trial{:}]};
    tmpData.time = {linspace(1 / tmpData.fsample, size(tmpData.trial{1}, 2) / tmpData.fsample, size(tmpData.trial{1}, 2))};
    tmpData = rmfield(tmpData, 'cfg');
    tmpData.grad.chanunit = repmat({'T'}, 1, size(tmpData.grad.chanpos, 1))';
    tmpData.grad.balance.reject.chanunitnew = repmat({'T'}, 1, size(tmpData.grad.chanpos, 1))';
    D = spm_eeg_ft2spm(tmpData, 'tmpftdata');
    
    S = [];
    S.D = D;
    D = spm_opm_amm(S);
    
    % Convert back to FieldTrip format
    tmpData = spm2fieldtrip(D);
    clear D;
    tmpData.grad.chantype = repmat({'MEG'}, size(tmpData.grad.chanpos, 1), 1);

    cfg = [];
    cfg.trl = [data_in.sampleinfo, zeros(size(data_in.sampleinfo, 1), 1)];
    data_in = ft_redefinetrial(cfg, tmpData);

    %% Step 9: Replace Bad Data with NaNs
    data_in = putbackNans(data_in, badGuidePost);
    opmSel = find(strncmp(data_out.label, 'G', 1));

    for tr = 1:numel(data_out.trial)
        data_out.trial{tr}(opmSel, :) = data_in.trial{tr};
    end

    % Fix the header to account for removed channels
    data_out = correctHeaderLabels(data_out);

    %% Step 10: Inspect Spectra After Processing (Optional Plotting)
    if plotop
        cfg = [];
        cfg.channel = data_in.label(opmSel);
        cfg.trial_length = 3;
        cfg.method = 'tim';
        cfg.foi = [1, 98];
        cfg.plot = 'yes';
        ft_opm_psd(cfg, data_in);
        title('After HFC');
        axis square;
        set(gcf, 'Position', [100.0000, -9.4000, 799.4000, 792.0000]);
        ylim([1e-1, 1e5]);
    end

end

%% Helper Functions
function data = correctHeaderLabels(data)
    cl_org = data.hdr.label;
    cl_new = data.label;
    [c, ia] = setdiff(cl_org, cl_new);
    data.hdr.label(ia) = [];
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

function comp = runICA(data_in)
    % Run ICA on the given data
    dat = cat(2, data_in.trial{:});
    dat(isnan(dat)) = 0;
    n_ic = fix(rank(dat));
    cfg = [];
    cfg.method = 'runica';
    cfg.numcomponent = n_ic;
    comp = ft_componentanalysis(cfg, data_in);
    comp.grad.chanunit = repmat({'T'}, 1, size(comp.grad.chanori, 1))';
    comp.grad.chantype = repmat({'MEG'}, 1, size(comp.grad.chanori, 1))';
    comp.grad.balance.comp.chanunitnew = repmat({'T'}, 1, size(comp.grad.balance.comp.chanunitnew, 1))';
    comp.grad.balance.comp.chantypenew = repmat({'MEG'}, 1, size(comp.grad.balance.comp.chanunitnew, 1))';
    
    % Visualize components
    cfg = [];
    cfg.laycfg.channel = '*Y';
    cfg.laycfg.grad = data_in.grad;
    cfg.layout = ft_prepare_layout(cfg.laycfg, data_in);
    cfg.viewmode = 'component';
    cfg.zlim = 'maxmin';
    cfg.compscale = 'local';
    cfg.contournum = 6;
    cfg.artifactalpha = 0.8;
    ft_databrowser(cfg, comp);
end