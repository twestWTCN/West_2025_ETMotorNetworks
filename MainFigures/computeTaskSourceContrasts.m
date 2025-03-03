function [] = computeTaskSourceContrasts(R, EEGSel, PATSel)
% computeTaskSourceContrasts: Compute task-related source contrasts for EEG or OPM data.
% This script runs a contrast analysis between different phases of task performance using DICS filters.

close all;
% Setup
cmap = brewermap(256, '*RdYlBu'); % Colormap setup
fresh = 0; % Set fresh to 1 if new computations are needed, otherwise load from saved data
pltflag = 1; % Flag to indicate whether to plot individual results
partit = {'Rest2Posture', 'Posture2Prep', 'Prep2Exec', 'Exec2Hold'}; % Task Phases

%% Choose Subjects
if EEGSel == 1
    R.import.type = 'EEG';
    R.import.site = 'DL';
    if PATSel == 1
        R.locSrcs.subsel = R.import.subselPat; % EEG patients
    elseif PATSel == 2
        R.locSrcs.subsel = R.import.subselCont; % EEG controls
    elseif PATSel == 3
        R.locSrcs.subsel = R.import.subselAll; % All EEG subjects
    end
else
    R.import.type = 'OPM';
    R.import.site = 'UCL';
    if PATSel == 1
        R.locSrcs.subsel = R.import.subselOPPat; % OPM patients
    elseif PATSel == 2
        R.locSrcs.subsel = R.import.subselOPCont; % OPM controls
    elseif PATSel == 3
        R.locSrcs.subsel = R.import.subselOPAll; % All OPM subjects
    end
end

%% Load Template Source Model and Template MRI
load(fullfile(R.path.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'), 'sourcemodel');
template_grid = sourcemodel;
clear sourcemodel;

% Load template MRI and atlas
templatefile = [R.path.ftpath '/template/anatomy/single_subj_T1.nii'];
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'spm';
atlas = ft_read_atlas([R.path.ftpath '/template/atlas/aal/ROI_MNI_V4.nii']);

% Define ROI locations
loclist = {[-37 -25 62], [-20 -75 44], [-16 -24 66], [35 -15 50]};

%% Define and Select Frequency Band
bnd = 2; % Set the desired frequency band index
bnddef = {[3 8], [14 30], [40 85], [2 98]};
bnddefname = {'3_8Hz', '14_30Hz', '40_85Hz'};

%% Main Loop Over Task Phases
for prt = 3 % Specify which task phase to process (e.g., prt=3 for Prep2Exec)
    subI = 0;
    for sub = R.locSrcs.subsel
        subI = subI + 1;
        if fresh
            %% Load Precomputed Filters and Leadfield
            headmodel = loadExpData(R, sub{1}, 'ForwardModel', '', 'headmodel', 'Mesh');
            sourceDICComFilt = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourceDICComFiltV3', 'inverse');
            leadfield = loadExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', 'solution');

            %% Retrieve Concatenated Data for Pre and Post Conditions
            partPair = {[2 2], [3 3], [4 4], [5 5]}; % Task part pairs for contrast analysis
            retdat = [];
            retdat.sub = sub;
            retdat.block = 3; % Task block
            retdat.cond = 1:4; % Conditions
            retdat.part = partPair{prt}(1); % Part before contrast
            retdat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']";
            retdat.subfold = "['Condition' num2str(cond) '_epoched']";
            retdat.ftflag = 1; % Flag for FT concatenation
            ftdata_C1 = retrieveData(R, retdat);

            retdat.part = partPair{prt}(2); % Part after contrast
            ftdata_C2 = retrieveData(R, retdat);

            % Remove short trials (< 512 samples)
            cfg = [];
            cfg.trials = cellfun(@length, ftdata_C1.trial) > 512;
            ftdata_C1 = ft_preprocessing(cfg, ftdata_C1);
            cfg.trials = cellfun(@length, ftdata_C2.trial) > 512;
            ftdata_C2 = ft_preprocessing(cfg, ftdata_C2);

            %% Define Pre and Post Epoch Details
            [pre, post] = getPrePostEpochDetails(prt + 1);

            % Extract Pre Data
            cfg = [];
            cfg.latency = pre;
            cfg.channel = intersect(ftdata_C1.label, ftdata_C2.label);
            ftdata_C1 = ft_selectdata(cfg, ftdata_C1);

            % Extract Post Data
            cfg.latency = post;
            ftdata_C2 = ft_selectdata(cfg, ftdata_C2);

            %% Remove NaNs from Trials
            nantr1 = ft_nancheck(ftdata_C1);
            nantr2 = ft_nancheck(ftdata_C2);
            cfg = [];
            cfg.trials = setdiff(1:numel(ftdata_C1.trial), nantr1);
            ftdata_C1 = ft_selectdata(cfg, ftdata_C1);
            cfg.trials = setdiff(1:numel(ftdata_C2.trial), nantr2);
            ftdata_C2 = ft_selectdata(cfg, ftdata_C2);

            %% Apply DICS Filter for Source Analysis
            comfilt = sourceDICComFilt{bnd}.avg;
            cfg = [];
            cfg.channel = intersect(ftdata_C1.label, comfilt.label);
            ftdata_C1 = ft_selectdata(cfg, ftdata_C1);
            ftdata_C2 = ft_selectdata(cfg, ftdata_C2);

            source_C1 = runBeamFormer(ftdata_C1, headmodel, leadfield, 'no', comfilt, 'DICS', bnddef{bnd}, R.import.type);
            source_C2 = runBeamFormer(ftdata_C2, headmodel, leadfield, 'no', comfilt, 'DICS', bnddef{bnd}, R.import.type);
            sourceDiff = source_C2;

            %% Calculate the Difference Between Post and Pre Conditions
            sourceDiff.avg.pow = (source_C2.avg.pow - source_C1.avg.pow) ./ ((source_C2.avg.pow + source_C1.avg.pow) / 2);

            %% Interpolate to Template MRI
            sourceDiff.pos = template_grid.pos;
            cfg = [];
            cfg.parameter = 'avg.pow';
            source_int = ft_sourceinterpolate(cfg, sourceDiff, template_mri);

            %% Flip Left/Right for Hand Dominance
            subTab = retrieveSubjectInfoTab(R, sub);
            source_int.pow = flipVolume(table2array(subTab(:, 'Hand')), source_int.pow, source_int.dim);

            %% Save Results as NIfTI
            statvol = [];
            statvol.anatomy = reshape(source_int.pow, source_int.dim);
            statvol.dim = source_int.dim;
            statvol.transform = source_int.transform;
            path = [R.path.datapath '\' R.path.expname '\Group\' partit{prt} '\PosturalHold\' sub{1} '_' partit{prt} '_' bnddefname{bnd}];
            ft_write_mri([path '.nii'], statvol, 'dataformat', 'nifti_spm');
            saveExpData(R, sub{1}, 'FunctionalVolumesV2', '', [partit{prt} '_sourcemapsDM_' bnddefname{bnd}], source_int, partit{prt});

            sourcesave{subI} = source_int;
        else
            source_int = loadExpData(R, sub{1}, 'FunctionalVolumesV2', '', [partit{prt} '_sourcemapsDM_' bnddefname{bnd}], partit{prt});
            sourcesave{subI} = source_int;
        end

        %% Plot Individual Results
        if pltflag
            h = figure(subI);
            cfg = [];
            cfg.method = 'ortho';
            cfg.funparameter = 'avg.pow';
            cfg.funcolorlim = [-0.5 0.5];
            cfg.atlas = atlas;
            cfg.figure = h;
            cfg.funcolormap = cmap;
            cfg.location = loclist{prt};
            ft_sourceplot(cfg, source_int);
            h.Name = sub{1};
        end
    end

    %% Group Statistics - Compute T-statistic
    for subI = 1:numel(sourcesave)
        sourcesave{subI}.pow = reshape(sourcesave{subI}.pow, sourcesave{subI}.dim);
        hsize = [10, 10, 10]; % Smoothing size
        sig = 1; % Smoothing sigma
        sourcesave{subI}.powSm = convolve3DGaussian(sourcesave{subI}.pow, hsize, sig);
        sourcesave{subI}.dim = [91, 109, 91];
        sourcesave{subI}.freq = 5;
    end

    stat = computeSourceTStat(sourcesave, 'powSm');
    df = numel(sourcesave) - 1;

    % Mask T-statistics
    stat.statmask = abs(stat.stat);
    stat.statmask(isnan(stat.stat)) = nan;

    % Plot T-statistic
    h = figure(200 + prt); clf;
    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'stat';
    cfg.funcolorlim = [-5 5];
    cfg.funcolormap = cmap;
    cfg.location = loclist{prt};
    cfg.atlas = atlas;
    cfg.figure = h;
    cfg.crosshair = 'off';
    ft_sourceplot(cfg, stat);
end
end