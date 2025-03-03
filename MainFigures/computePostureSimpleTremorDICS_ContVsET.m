function [] = computePostureSimpleTremorDICS_ContVsET(R, EEGSel)
% computePostureSimpleTremorDICS_ContVsET: Compare synchronized tremor sources ET patients and control subjects
% using coherence analysis in the postural hold condition. This function computes DICS (Dynamic Imaging of Coherent Sources)
% for EEG or OPM data and compares ET versus control groups.

if EEGSel == 1
    R.locSrcs.subsel = R.import.subselAll; % EEG
else
    R.locSrcs.subsel = R.import.subselOPAll; % OPM
end

close all;
refch = 'pca001'; % Reference channel for computing coherence

subI = 0;

% Load Template Grid
load(fullfile(R.path.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel;

templatefile = [R.path.ftpath '/template/anatomy/single_subj_T1.nii'];
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'spm';

atlas = ft_read_atlas([R.path.ftpath '/template/atlas/aal/ROI_MNI_V4.nii']);

% Colormap settings
cmap = brewermap(512, '*RdYlBu');
cmap = cmap(257:512, :);
loclist = {[-18 -6 68]}; % Specify location for further analyses

fresh = 1;
for sub = R.locSrcs.subsel
    subI = subI + 1;
    if fresh
        block = 2; % Posture task
        
        % Load precomputed filters and leadfield
        headmodel = loadExpData(R, sub{1}, 'ForwardModel', '', 'headmodel', 'Mesh');
        sourceDICComFilt = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourceDICRefACCComFiltV3', 'inverse');
        leadfield = loadExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', 'solution');

        % Retrieve postural data
        retdat = [];
        retdat.sub = sub;
        retdat.block = block;
        retdat.fileapend = "'pp_arV3'";
        retdat.subfold = 'epoched';
        retdat.ftflag = 1; % FT concatenation flag
        ftdata_posture = retrieveData(R, retdat);

        % Filter trials with excessive movement
        cfg = [];
        cfg.channel = {'AccX', 'AccY', 'AccZ'};
        accData = ft_selectdata(cfg, ftdata_posture);
        stdtr = cellfun(@(x) std(x')', accData.trial, 'UniformOutput', 0);
        sel = find(abs(zscore(max([stdtr{:}]))) < 3);
        cfg = [];
        cfg.trials = sel;
        ftdata_posture = ft_selectdata(cfg, ftdata_posture);

        % Get peak tremor frequency from PCA
        [~, peakFrq, ~, ~, ~, accData] = computeTremorProperties(ftdata_posture, [], 0);
        cfg = [];
        ftdata_posture = ft_appenddata(cfg, ftdata_posture, accData);

        % Intersect channels
        cfg = [];
        cfg.channel = sourceDICComFilt{1}.avg.label;
        cfg.channel = [cfg.channel; refch];
        ftdata_posture = ft_selectdata(cfg, ftdata_posture);

        % Demean the data
        cfg = [];
        cfg.demean = 'yes';
        ftdata_posture = ft_preprocessing(cfg, ftdata_posture);

        % Apply the DICS filter
        bnddef{1} = [peakFrq - 1.5, peakFrq + 1.5];
        bnd = 1;
        comfilt = sourceDICComFilt{bnd}.avg;
        source_POST = runBeamFormer(ftdata_posture, headmodel, leadfield, 'no', comfilt, 'DICSREF', bnddef{bnd}, R.import.type, {refch});
        source_POST.pos = template_grid.pos;

        % Interpolate with template MRI
        cfg = [];
        cfg.parameter = 'avg.coh';
        source_int = ft_sourceinterpolate(cfg, source_POST, template_mri);

        % Z-normalize coherence
        edf = 1 / (size(source_POST.cumtapcnt, 1) - 2); % Effective degrees of freedom
        source_int.cohz = (atanh(source_int.coh) - edf) / sqrt(edf);

        % Flip LR if required (for left/right-handers)
        subTab = retrieveSubjectInfoTab(R, sub);
        source_int.coh = flipVolume(table2array(subTab(:, 'Hand')), source_int.coh, source_int.dim);
        source_int.cohz = flipVolume(table2array(subTab(:, 'Hand')), source_int.cohz, source_int.dim);

        % Convert atlas to the correct coordinate system
        atlas = ft_convert_coordsys(atlas, source_int.coordsys);

        % Plot source results
        h = figure(subI); clf;
        cfg = [];
        cfg.method = 'ortho';
        cfg.funparameter = 'cohz';
        cfg.funcolorlim = [-0.1, 0.1];
        cfg.location = [-46, -20, 64];
        cfg.funcolormap = cmap;
        cfg.atlas = atlas;
        cfg.figure = h;
        ft_sourceplot(cfg, source_int);

        % Find local source peaks
        roi_center = [-6, 0, 70];
        radius = 1000; % mm for global peak
        global_peak = findLocalSourcePeak(source_int, roi_center, radius, 'cohz');
        radius = 50; % mm for local peak
        local_peak = findLocalSourcePeak(source_int, roi_center, radius, 'cohz');
        maxLabel = atlas_lookup_PublicFT(atlas, local_peak.pos, 'coordsys', 'spm');
        locSave(:, subI) = {local_peak.pos, local_peak.value, global_peak.pos, global_peak.value, peakFrq, maxLabel{1}};

        % Save the volume data
        pathsave = saveExpData(R, sub{1}, 'DICSImages', '', [sub{1}, 'PosturalHold_tremorDICS_'], source_int, 'Tremor');
        sourcesave{subI} = source_int;
    else
        sourcesave{subI} = loadExpData(R, sub{1}, 'DICSImages', '', [sub{1}, 'PosturalHold_tremorDICS_'], 'Tremor');
    end
end

% Process saved data
for subI = 1:numel(sourcesave)
    sourcesave{subI}.coh = reshape(sourcesave{subI}.coh, sourcesave{subI}.dim);
    sourcesave{subI}.cohz = reshape(sourcesave{subI}.cohz, sourcesave{subI}.dim);
    sourcesave{subI}.dim = [91, 109, 91];
    sourcesave{subI}.freq = 5;
end

% Split data into Control and ET groups
if EEGSel == 1
    % EEG List
    sourceCONT = sourcesave(1:11);
    sourceET = sourcesave([12:14, 16:23]); % Subj 15 is missing ACC
else
    % OPM List
    sourceCONT = sourcesave([1:4, 8]);
    sourceET = sourcesave([5:7, 9]);
end

% Compute Grand Averages for Control and ET groups
cfg = [];
cfg.parameter = 'cohz';
sourceGACONT = ft_sourcegrandaverage(cfg, sourceCONT{:});
sourceGAET = ft_sourcegrandaverage(cfg, sourceET{:});

% Calculate Difference between ET and Control Groups
sourceDIFF = sourcesave{1};
sourceDIFF.cohz = (sourceGAET.cohz - sourceGACONT.cohz);
sourceDIFF.coordsys = sourcesave{1}.coordsys;
sourceDIFF.anatomy = sourcesave{1}.anatomy;

% Plot the Difference
h = figure(100 + 1); clf;
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'cohz';
cfg.funcolormap = cmap;
cfg.atlas = atlas;
cfg.figure = h;
ft_sourceplot(cfg, sourceDIFF);

%% Statistical Analysis - Compute T-statistics
% Calculate T-statistics between ET and Control Groups
stat = computeSourceTStatPaired(sourceET, sourceCONT, 'cohz');
stat.stat(abs(stat.stat) < tinv(0.95, numel(sourcesave) - 2)) = 0;
stat.statmask = abs(stat.stat) > tinv(0.95, numel(sourcesave) - 2);
stat.statmask(isnan(stat.statmask)) = 0;

% Plot T-statistics
h = figure(200); clf;
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'stat';
cfg.maskparameter = 'stat';
cfg.maskstyle = 'opacity';
cfg.opacitylim = [0, 1];
cfg.funcolorlim = [-3, 3];
cfg.funcolormap = cmap;
cfg.location = loclist{1};
cfg.atlas = atlas;
cfg.figure = h;
cfg.crosshair = 'no';
cfg.axis = 'off';
ft_sourceplot(cfg, stat);
end
