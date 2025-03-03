function prepareMRI_OPM_normalizedGrid(R)
% Prepare MRI and forward models for OPM data Source analysis
% Input: R - structure with necessary paths and subject details

close all;
subI = 0;

for sub = R.locSrcs.subsel
    subI = subI + 1;
    flag_check = true;
    interactiveRA = false;
    fresh = true;
    
    if fresh
        freshflag = true;
    else
        try
            disp('Trying to load headmodel');
            headmodel = loadExpData(R, sub{1}, 'ForwardModel', '', 'headmodel', 'Mesh');
            disp('Trying to load leadfield');
            leadfield = loadExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', 'solution');
            freshflag = false;
        catch
            freshflag = true;
        end
    end
    
    if freshflag
        try
            disp('Loading pre-segmented MRI');
            mri_segmented = loadExpData(R, sub{1}, 'ForwardModel', '', 'segmented', 'mri');
            mri = loadExpData(R, sub{1}, 'ForwardModel', '', 'importmri', 'mri');
        catch
            disp('Couldn''t find segmented MRI, segmenting now. This may take a few minutes...');
            mri = ft_read_mri(fullfile(R.path.datapath, R.path.expname, sub{1}, 'MRI', ['wanon_' sub{1} '_t1_mpr.nii']), 'dataformat', 'nifti_spm');
            mri.coordsys = 'mni';
            saveExpData(R, sub{1}, 'ForwardModel', '', 'importmri', mri, 'mri');
            
            cfg = [];
            cfg.write = 'no';
            cfg.spmversion = 'spm12';
            cfg.output = 'brain';
            mri_segmented = ft_volumesegment(cfg, mri);
            saveExpData(R, sub{1}, 'ForwardModel', '', 'segmented', mri_segmented, 'mri');
        end
        
        cfg = [];
        cfg.spmversion = 'spm12';
        cfg.method = 'singleshell';
        cfg.tissue = 'brain';
        headmodel = ft_prepare_headmodel(cfg, mri_segmented);
        
        dummyft = loadExpData(R, sub{1}, 'Rest', [], 'rawdata', 'raw'); % Load dummy data for gradient structure
        
        if interactiveRA
            grad = dummyft.grad;
            elec = [];
            elec.elecpos = grad.coilpos;
            elec.chanpos = grad.coilpos;
            elec.label = grad.label;
            
            cfg = [];
            cfg.method = 'interactive';
            cfg.headshape = headmodel.bnd;
            gradRealigned = ft_electroderealign(cfg, elec);
            gradRealigned2 = ft_transform_geometry(gradRealigned.homogeneous, grad);
        else
            gradRealigned2 = dummyft.grad;
        end
        
        % Plot sensor and headmodel
        figure;
        ft_plot_sens(gradRealigned2);
        hold on;
        ft_plot_headmodel(headmodel);
        
        % Construct template matched grid
        load(fullfile(R.path.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'));
        load(fullfile(R.path.ftpath, 'template/headmodel/standard_bem'));
        template_grid = sourcemodel;
        template_grid = ft_convert_units(template_grid, 'mm');
        template_grid.inside = ft_inside_headmodel(template_grid.pos, vol, 'inwardshift', 0);
        clear sourcemodel vol;
        
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
        
        % Plot the subject headmodel and grid positions
        figure; hold on;
        ft_plot_headmodel(headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
        ft_plot_mesh(grid.pos(grid.inside, :));
        
        cfg = [];
        cfg.grad = gradRealigned2;
        cfg.headmodel = headmodel;
        cfg.grid = grid;
        cfg.channel = dummyft.label(strncmp(dummyft.label, 'G', 1));
        cfg.grid.unit = 'mm';
        cfg.normalize = 'yes';
        leadfield = ft_prepare_leadfield(cfg);
        
        if flag_check
            sanity_check_leadfield(gradRealigned2, leadfield, headmodel, headmodel.bnd(1));
        end
        
        % Save computed data
        saveExpData(R, sub{1}, 'ForwardModel', '', 'headmodel', headmodel, 'Mesh');
        saveExpData(R, sub{1}, 'ForwardModel', '', 'sourcemodel', grid, 'solution');
        saveExpData(R, sub{1}, 'ForwardModel', '', 'leadfieldNew', leadfield, 'solution');
    end % End Loading/Computing Source Model
end

end