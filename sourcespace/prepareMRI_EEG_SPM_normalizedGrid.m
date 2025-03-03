function prepareMRI_EEG_SPM_normalizedGrid(R)
    % Function to prepare MRI and EEG data with SPM and Fieldtrip for a normalized grid
    % Arguments:
    %   R - Structure containing paths and configurations
    %
    % This function prepares head models, sourcemodels, realigns electrodes,
    % and creates leadfields for subjects selected in the 'R' structure.

    % Initialization
    fresh = 0;
    templateMRI = fullfile(R.path.ftpath, 'template', 'anatomy', 'single_subj_T1.nii');
    subI = 0;
    
    % Iterate over each subject in the selected list
    for sub = R.locSrcs.subsel
        subI = subI + 1;
        flag_check = 1;
        freshflag = determine_freshflag(R, sub{1}, fresh);
        
        %% Construct HeadModel using SPM
        if freshflag
            [headmodel, elec_realigned] = construct_headmodel(R, sub{1}, templateMRI);
        else
            headmodel = loadExpData(R, sub{1}, 'ForwardModel', '', 'headmodel', 'Mesh');
        end
        
        if flag_check
            visualize_headmodel(headmodel);
        end

        %% Create Sourcemodel for the subject
        freshflag = determine_freshflag(R, sub{1}, fresh);
        if freshflag
            grid = create_sourcemodel(R, sub{1}, templateMRI);
        else
            grid = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourcemodel', 'solution');
        end
        
        if flag_check
            visualize_sourcemodel(headmodel, grid);
        end

        %% Realign Electrodes
        freshflag = determine_freshflag(R, sub{1}, fresh);
        if freshflag
            elec_realigned = realign_electrodes(R, sub{1}, headmodel);
        else
            elec_realigned = loadExpData(R, sub{1}, 'ForwardModel', '', 'elec_realigned', 'solution');
        end
        
        if flag_check
            visualize_electrodes(headmodel, elec_realigned);
        end

        %% Create Leadfield
        freshflag = determine_freshflag(R, sub{1}, fresh);
        if freshflag
            leadfield = create_leadfield(headmodel, elec_realigned, grid);
        else
            leadfield = loadExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', 'solution');
        end
        
        %% Check Results
        if flag_check
            sanity_check_leadfield(elec_realigned, leadfield, headmodel, headmodel.bnd(1));
        end
        
        disp(['Subject ', sub{1}, ' forward model complete']);
        close all;
    end
end

%% Helper Functions

function freshflag = determine_freshflag(R, subject, fresh)
    % Determine if fresh computation is needed
    if fresh
        freshflag = 1;
    else
        try
            loadExpData(R, subject, 'ForwardModel', '', 'headmodel', 'Mesh');
            freshflag = 0;
        catch
            freshflag = 1;
        end
    end
end

function [headmodel, elec_realigned] = construct_headmodel(R, subject, templateMRI)
    % Load electrode definitions and prepare SPM headmodel
    elec_default = ft_read_sens(fullfile(R.path.ftpath, 'template', 'electrode', 'standard_1005.elc'));
    elec_default = ft_convert_units(elec_default, 'mm');
    
    fttmp = loadExpData(R, subject, 'Rest', [], 'rawdata', 'raw');
    cfg = [];
    cfg.trials = cellfun(@length, fttmp.trial) > 512;
    fttmp = ft_preprocessing(cfg, fttmp);
    fttmp = makeSameTime(fttmp);
    
    cfg = [];
    cfg.latency = [0 1];
    fttmp = ft_selectdata(cfg, fttmp);
    
    cfg = [];
    cfg.channel = intersect(fttmp.label, elec_default.label);
    fttmp = ft_selectdata(cfg, fttmp);
    fttmp = makeSameTime(fttmp);
    
    spmDataPath = fullfile(R.path.datapath, R.path.expname, subject, 'SPM', 'SPMdummy.mat');
    mkdir(fileparts(spmDataPath));
    D = spm_eeg_ft2spm(fttmp, spmDataPath);
    
    D = sensors(D, 'EEG', elec_default);
    L.fid.pnt = elec_default.chanpos(1:3, :);
    L.fid.label = {'LPA', 'RPA', 'Nasion'};
    D = fiducials(D, L);
    save(D);
    clear fttmp;
    
    MRIPath = fullfile(R.path.datapath, R.path.expname, subject, 'MRI', ['anon_', subject, '_t1_mpr.nii']);
    if ~exist(MRIPath, 'file')
        MRIPath = [];
        FidTab = MRIFiducials;
        Fid = FidTab(find(contains(FidTab(:, 1), 'TEMPLATE')), 2:end);
    else
        FidTab = MRIFiducials;
        Fid = FidTab(find(contains(FidTab(:, 1), subject)), 2:end);
    end
    makeSPMHeadModel(R, spmDataPath, MRIPath, Fid);
    
    D = spm_eeg_load(spmDataPath);
    A = spm_eeg_inv_get_vol_sens(D);
    clear D;
    
    headmodel = A.EEG.vol;
    headmodel = ft_convert_units(headmodel, 'mm');
    saveExpData(R, subject, 'ForwardModel', '', 'headmodel', headmodel, 'Mesh');
end

function visualize_headmodel(headmodel)
    figure; hold on;
    ft_plot_headmodel(headmodel, 'facealpha', 0.5);
    view(90, 0);
end

function grid = create_sourcemodel(R, subject, templateMRI)
    try
        MRIPath = fullfile(R.path.datapath, R.path.expname, subject, 'MRI', ['wanon_', subject, '_t1_mpr.nii']);
        mri = ft_read_mri(MRIPath);
    catch
        mri = ft_read_mri(templateMRI);
    end
    
    load(fullfile(R.path.ftpath, 'template', 'sourcemodel', 'standard_sourcemodel3d10mm'));
    load(fullfile(R.path.ftpath, 'template', 'headmodel', 'standard_bem'));
    
    template_grid = sourcemodel;
    template_grid = ft_convert_units(template_grid, 'mm');
    template_grid.inside = ft_inside_headmodel(template_grid.pos, vol, 'inwardshift', 0);
    clear vol sourcemodel;
    
    mri.coordsys = 'ras';
    cfg = [];
    cfg.warpmni = 'yes';
    cfg.template = template_grid;
    cfg.nonlinear = 'yes';
    cfg.mri = mri;
    cfg.unit = 'mm';
    cfg.spmversion = 'spm12';
    cfg.spmmethod = 'new';
    grid = ft_prepare_sourcemodel(cfg);
    
    saveExpData(R, subject, 'ForwardModel', '', 'sourcemodel', grid, 'solution');
end

function visualize_sourcemodel(headmodel, grid)
    figure; hold on;
    ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
    ft_plot_mesh(grid.pos(grid.inside, :));
    view(125, 10);
end

function elec_realigned = realign_electrodes(R, subject, headmodel)
    elec_default = ft_read_sens(fullfile(R.path.ftpath, 'template', 'electrode', 'standard_1005.elc'));
    elec_default = ft_convert_units(elec_default, 'mm');
    
    scalp_index = determine_scalp_index(headmodel);
    cfg = [];
    cfg.method = 'project';
    cfg.headshape = headmodel.bnd(scalp_index);
    elec_realigned = ft_electroderealign(cfg, elec_default);
    elec_realigned = ft_convert_units(elec_realigned, 'mm');
    
    saveExpData(R, subject, 'ForwardModel', '', 'elec_realigned', elec_realigned, 'solution');
end

function visualize_electrodes(headmodel, elec_realigned)
    figure; hold on;
    ft_plot_sens(elec_realigned, 'elecsize', 40);
    ft_plot_headmodel(headmodel, 'facealpha', 0.5);
    view(90, 0);
end

function scalp_index = determine_scalp_index(headmodel)
    if strcmp(headmodel.type, 'bemcp')
        scalp_index = 3;
    else
        scalp_index = 1;
    end
end

function leadfield = create_leadfield(headmodel, elec_realigned, grid)
    cfg = [];
    cfg.elec = elec_realigned;
    cfg.headmodel = headmodel;
    cfg.grid = grid;
    cfg.sourcemodel.unit = 'mm';
    cfg.normalize = 'yes';
    leadfield = ft_prepare_leadfield(cfg);
    saveExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', leadfield, 'solution');
end