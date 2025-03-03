function [] = computePostureSimpleTremorDICS_Overlap(R, EEGSel)
% computePostureSimpleTremorDICS_Overlap: Compute overlap of tremor DICS images for postural hold condition.
% This function computes coherence-based tremor analysis using DICS (Dynamic Imaging of Coherent Sources).
% It compares EEG versus OPM data for subjects in a postural hold condition and produces overlap maps.

if EEGSel == 1
    R.locSrcs.subsel = R.import.subselPat; % EEG
else
    R.locSrcs.subsel = R.import.subselOPPat; % OPM
end

close all;
refch = 'pca001'; % Reference channel used for computing coherence

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

% Specify location for further analyses
loclist = {[-18 -6 68]};

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

        % Save the volume data
        saveExpData(R, sub{1}, 'DICSImages', '', [sub{1}, 'PosturalHold_tremorDICS_'], source_int, 'Tremor');
        sourcesave{subI} = source_int;
    else
        sourcesave{subI} = loadExpData(R, sub{1}, 'DICSImages', '', [sub{1}, 'PosturalHold_tremorDICS_'], 'Tremor');
    end
end

% Process saved data
for subI = 1:numel(sourcesave)
    coh = sourcesave{subI}.cohz;
    ZCoh = coh;
    sourcesave{subI}.coh = reshape(sourcesave{subI}.coh, sourcesave{subI}.dim);
    sourcesave{subI}.cohz = reshape(sourcesave{subI}.cohz, sourcesave{subI}.dim);
    sourcesave{subI}.ZCOH = reshape(ZCoh, sourcesave{subI}.dim);
    hsize = [10, 10, 10];
    sig = 3;
    sourcesave{subI}.ZCOHsm = convolve3DGaussian(sourcesave{subI}.ZCOH, hsize, sig);

    sourcesave{subI}.dim = [91, 109, 91];
    sourcesave{subI}.freq = 5;
end

% Subject selection based on EEG or OPM
if EEGSel == 1
    % EEG List
    leaveOut = [4]; % Subject 4 doesnt have accelerometer
    list = setdiff(1:12, leaveOut);
else
    % OPM List
    leaveOut = [];
    list = setdiff(1:4, leaveOut);
end

%% Compute Overlap Figure
olap = [];
for subI = 1:numel(list)
    olap(:, :, :, subI) = sourcesave{list(subI)}.cohz > prctile(sourcesave{list(subI)}.cohz(:), 85);
end

olap = sum(olap, 4);
threshmap = (olap > 1) & (sourcesave{subI}.inside == 1);

overlapSource = sourcesave{1};
overlapSource.olap = olap;
overlapSource.threshmap = olap / max(olap(:));
cmapSet = brewermap(fix(numel(list) / 1), 'RdBu');
h = figure(400 + 1); clf;
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'olap';
cfg.funcolormap = cmapSet;
cfg.funcolorlim = [0, numel(list)];
cfg.maskparameter = 'threshmap';
cfg.maskstyle = 'opacity';
cfg.atlas = atlas;
cfg.figure = h;
cfg.crosshair = 'off';
cfg.location = loclist{1};
ft_sourceplot(cfg, overlapSource);
end
