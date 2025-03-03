function localizeSourcesCommonFilter_normalizedGrid(R)
% Localize sources using a common filter and normalized grid
% Input: R - structure with necessary paths and subject details

close all;

% Load Template Grid
load(fullfile(R.path.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel;

% Load Template MRI
templatefile = fullfile(R.path.ftpath, 'template/anatomy/single_subj_T1.nii');
template_mri = ft_read_mri(templatefile);

subI = 0;
for sub = R.locSrcs.subsel
    subI = subI + 1;
    fresh = true;
    
    if fresh
        freshflag = true;
    else
        try
            sourceDICComFilt = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourceDICComFilt', 'inverse');
            sourceLCMVComFilt = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourceLCMVComFilt', 'inverse');
            freshflag = false;
        catch
            freshflag = true;
        end
    end

    if freshflag
        % Retrieve task data
        retdat = struct();
        retdat.sub = sub;
        retdat.block = 3; % task block
        retdat.cond = 1:4; % condition
        retdat.part = 2:5; % movement parts
        retdat.fileapend = "[R.epoch.names{part} '_block_pp_arV3']";
        retdat.subfold = "['Condition' num2str(cond) '_epoched']";
        retdat.ftflag = 1; % Use FieldTrip concatenation
        ftdata_cat_TSK = retrieveData(R, retdat);

        retdat.block = 1:2; % rest blocks
        retdat.fileapend = "'pp_arV3'";
        retdat.subfold = 'epoched';
        retdat.ftflag = 1;
        ftdata_cat_RP = retrieveData(R, retdat);

        % Ensure trial info is matched
        maxTr = 6;
        if isfield(ftdata_cat_RP, 'trialinfo')
            ftdata_cat_RP.trialinfo = fillTrialInfo(ftdata_cat_RP.trialinfo, maxTr);
        end

        % Append task and rest data
        cfg = [];
        cfg.keepsampleinfo = 'no';
        cfg.appenddim = 'rpt';
        ftdataCat = ft_appenddata(cfg, ftdata_cat_TSK, ftdata_cat_RP);

        % Preprocess and ensure uniform time length
        cfg = [];
        cfg.trials = cellfun(@length, ftdataCat.trial) > 1024;
        ftdataCat = ft_preprocessing(cfg, ftdataCat);
        ftdataCat = makeSameTime(ftdataCat);

        cfg.latency = [0 2];
        ftdataCat = ft_selectdata(cfg, ftdataCat);
        ftdataCat = makeSameLength(ftdataCat, 1024);

        % Load precomputed headmodel and leadfield
        try
            disp('Loading headmodel and leadfield');
            headmodel = loadExpData(R, sub{1}, 'ForwardModel', '', 'headmodel', 'Mesh');
            leadfield = loadExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', 'solution');
        catch
            error(['No precomputed forward model found for ' sub{1}]);
        end

        % Check coregistration
        figure;
        try
            elec_realigned = loadExpData(R, sub{1}, 'ForwardModel', '', 'elec_realigned', 'solution');
            ft_plot_sens(elec_realigned, 'elecsize', 50);
        catch
            ft_plot_sens(ftdataCat.grad);
        end

        ft_nancheck(ftdataCat);
        ftdataCat = getFTAccPca(ftdataCat, [3 8]);

        % Run beamforming for different frequency bands
        bnddef = {[3 8], [14 30], [40 85], [2 98]};
        sourceDICComFilt = cell(1, 4);
        sourceDICREFComFilt = cell(1, 4);
        for bnd = 1:4
            sourceDICComFilt{bnd} = runBeamFormer(ftdataCat, headmodel, leadfield, 'no', [], 'DICS', bnddef{bnd}, R.import.type, [], []);
            sourceDICREFComFilt{bnd} = runBeamFormer(ftdataCat, headmodel, leadfield, 'no', [], 'DICSREF', bnddef{bnd}, R.import.type, {'pca001'}, []);
            sourceDICComFilt{bnd}.pos = template_grid.pos;
            sourceDICREFComFilt{bnd}.pos = template_grid.pos;
        end
        saveExpData(R, sub{1}, 'ForwardModel', '', 'sourceDICComFiltV3', sourceDICComFilt, 'inverse');
        saveExpData(R, sub{1}, 'ForwardModel', '', 'sourceDICRefACCComFiltV3', sourceDICREFComFilt, 'inverse');

        % Run LCMV beamformer
        sourceLCMVComFilt = runBeamFormer(ftdataCat, headmodel, leadfield, 'no', [], 'LCMV', [], R.import.type, [], []);
        sourceLCMVComFilt.pos = template_grid.pos;
        saveExpData(R, sub{1}, 'ForwardModel', '', 'sourceLCMVComFiltV3', sourceLCMVComFilt, 'inverse');
    end

    % Load Example to Plot
    sourceDICComFilt = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourceDICComFiltV3', 'inverse');

    % The rest is intended for plotting (omitted)
end

end
