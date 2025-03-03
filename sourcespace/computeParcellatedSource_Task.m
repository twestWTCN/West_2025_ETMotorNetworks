function computeParcellatedSource_Task(R)
% Compute parcellated source for the task-related data
% Input: R - structure containing paths and subject details

close all;

% Load Template Grid
load(fullfile(R.path.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'));
template_grid = sourcemodel;

subI = 0;
for sub = R.locSrcs.subsel
    subI = subI + 1;
    for part = 2:5
        % Load precomputed filters
        sourceLCMVComFilt = loadExpData(R, sub{1}, 'ForwardModel', '', 'sourceLCMVComFiltV3', 'inverse');

        for cond = 1:4
            % Retrieve concatenated data
            retdat = struct();
            retdat.sub = sub;
            retdat.block = 3; % task block
            retdat.cond = cond; % condition
            retdat.part = part; % movement (prep, exec, hold)
            retdat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']";
            retdat.subfold = "['Condition' num2str(cond) '_epoched']";
            retdat.ftflag = 1; % Use FieldTrip concatenation
            ftdataCat = retrieveData(R, retdat);

            trlInfo = ftdataCat.trialinfo;
            % Remove trial info from data
            if isfield(ftdataCat, 'trialinfo')
                ftdataCat = rmfield(ftdataCat, 'trialinfo');
            end

            comfilt = sourceLCMVComFilt.avg;

            % Ensure data matches filter labels
            cfg = [];
            cfg.channel = intersect(ftdataCat.label, comfilt.label);
            ftdata = ft_selectdata(cfg, ftdataCat);

            source = [];
            source.pos = template_grid.pos;

            % Load atlas
            atlas = ft_read_atlas(fullfile(R.path.ftpath, 'template/atlas/aal/ROI_MNI_V4.nii'));
            
            % Interpolate atlas onto source model
            cfg = [];
            cfg.interpmethod = 'nearest';
            cfg.parameter = 'tissue';
            atlas2 = ft_sourceinterpolate(cfg, atlas, sourcemodel);

            % Get the inside positions of the atlas
            atlas2.inside = reshape(sourcemodel.inside, atlas2.dim);

            % Prepare virtual channels data structure
            ftdata_VC = [];
            ftdata_VC.fsample = ftdata.fsample;
            ftdata.time = ftdata.time;

            atlastis = atlas2.tissue(:);
            dS = 5; % Resolution
            c = 0; labellist = {}; newPos = []; newFilter = {};
            for i = 1:dS:numel(atlastis)
                if atlas2.inside(i) == 1
                    if atlastis(i) > 0 && ~isempty(comfilt.filter{i})
                        c = c + 1;
                        labellist{c} = [atlas2.tissuelabel{atlastis(i)} '_' num2str(c)];
                        newPos(c, :) = source.pos(i, :);
                        newFilter{c} = comfilt.filter{i};
                    end
                end
            end

            % Construct virtual channels
            VCData = constructVirtualChannels(newFilter, ftdata, 1:numel(ftdata.label), labellist);
            
            % Select non-sensor data based on modality type
            if strcmp(R.import.type, 'EEG')
                assert(~checkHeader(ftdataCat), 'Your header is corrupt!');
                nonSensInd = find(~contains(ftdataCat.hdr.chantype, 'EEG'));
            elseif strcmp(R.import.type, 'OPM')
                nonSensInd = find(~contains(ftdataCat.label, 'G'));
            end
            cfg = [];
            cfg.channel = ftdataCat.label(nonSensInd);
            nonSensData = ft_selectdata(cfg, ftdataCat);
            if isfield(nonSensData, 'elec'); nonSensData = rmfield(nonSensData, 'elec'); end
            if isfield(nonSensData, 'hdr');  nonSensData = rmfield(nonSensData, 'hdr'); end
            if isfield(nonSensData, 'grad');  nonSensData = rmfield(nonSensData, 'grad'); end

            % Concatenate virtual channels with non-sensor data
            VCDataCat = ft_appenddata([], VCData, nonSensData);
            VCDataCat.pos = [newPos; nan(size(VCDataCat.label, 1) - size(newPos, 1), 3)];
            VCDataCat.hdr.chantype = [repmat({'vc'}, 1, size(newPos, 1)), ftdataCat.hdr.chantype(nonSensInd)'];
            VCDataCat.hdr.label = VCDataCat.label;
            VCDataCat.hdr.nChans = numel(VCDataCat.label);
            VCDataCat.hdr.chanunit = nan(1, numel(VCDataCat.label));

            % Put back trial info
            VCDataCat.trialinfo = trlInfo;

            % Check consistency
            checkConsistentDataLength(VCDataCat, 0.1);
            assert(~checkHeader(VCDataCat), 'Your header is corrupt!');
            
            % Save processed data
            fileappend = [R.epoch.names{part} '_trans_pp_VC'];
            saveExpData(R, sub{1}, 'Task', [], fileappend, VCDataCat, ['Condition' num2str(cond) '_epoched']);
        end
    end
end
end
