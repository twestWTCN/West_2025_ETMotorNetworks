function [] = computePostureSimpleTremorDICS_RegressVC(R,EEGSel)
if EEGSel == 1
    R.locSrcs.subsel = R.import.subselAll; % EEG
else
    R.locSrcs.subsel = R.import.subselOPAll; % OPM
end

close all
% R.import.blockpart = R.locSrcs.blockpart;

% refch = 'EMG';
refch = 'pca001';
% refch = 'AccMag';

subI = 0;
% Load Template Grid
load(fullfile(R.path.ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm'));
template_grid = sourcemodel;
clear sourcemodel

templatefile = [R.path.ftpath '/template/anatomy/single_subj_T1.nii'];
template_mri = ft_read_mri(templatefile);
template_mri.coordsys = 'spm';

atlas = ft_read_atlas([R.path.ftpath '/template/atlas/aal/ROI_MNI_V4.nii']);

cmap = brewermap(512,'*RdYlBu');
cmap = cmap(257:512,:);

fresh = 0;
for sub =R.locSrcs.subsel
    subI = subI + 1;
    if fresh
        block = 2; % posture
        % Load precomputed filters
        headmodel = loadExpData(R,sub{1},'ForwardModel','','headmodel','Mesh');
        sourceDICComFilt = loadExpData(R,sub{1},'ForwardModel','','sourceDICRefACCComFiltV3','inverse');
        leadfield = loadExpData(R,sub{1},'ForwardModel','','leadfieldNew','solution');

        % Concatanated data
        retdat = [];
        retdat.sub = sub;
        retdat.block = block; % task
        retdat.fileapend = "'pp_arV3'";
        retdat.subfold = 'epoched';
        retdat.ftflag = 1; % FT concat
        ftdata_posture = retrieveData(R,retdat);

        %% Find trial by trial coherence and take out top 50%
        [~, peakFrq,~,~,~,accData] = computeTremorProperties(ftdata_posture,[],0);

        %% getReference Signal to regress out
        % Concatanated data
        retdat = [];
        retdat.sub = sub;
        retdat.block = block; % task
        retdat.fileapend = "'_pp_VC'";
        retdat.subfold = 'epoched';
        retdat.ftflag = 1; % FT concat
        vc_data = retrieveData(R,retdat);
        vc_data = rmfield(vc_data,'pos');

        ref_data = vc_data;
        subTab = retrieveSubjectInfoTab(R,sub);

        cfg = [];
        if strcmp(subTab.Hand{1},'L')
            cfg.channel = {'Supp_Motor_Area_R*','Postcentral_R*','Precentral_R*','Paracentral_Lobule_R*'};
        else
            cfg.channel = {'Supp_Motor_Area_L*','Postcentral_L*','Precentral_L*','Paracentral_Lobule_L*'};
        end
        cfg.method       = 'pca';
        ref_data = ft_componentanalysis(cfg,ref_data);
        ref_data = rmfield(ref_data,'topo');
        ref_data = rmfield(ref_data,'topolabel');
        ref_data = rmfield(ref_data,'unmixing');

        cfg = [];
        cfg.channel = 'pca001';
        Y1 = ft_selectdata(cfg,ref_data);
        cfg.channel = 'pca002';
        Y2 = ft_selectdata(cfg,ref_data);
        cfg.channel = 'pca003';
        Y3 = ft_selectdata(cfg,ref_data);

        % Select on coherence level.
        cx = []; cmag = [];
        for tr = 1:numel(accData.trial)
            [cx(:,tr),fx] = mscohere(accData.trial{tr},Y1.trial{tr},256,128,256,512);
            fxind = find(fx>=4 & fx<=12);
            cmag(tr) = max(cx(fxind,tr));
        end
        trsel = find(cmag>prctile(cmag,50));

        cfg.channel = {'eeg','O9','O10'};
        X = ft_selectdata(cfg,ftdata_posture);
        cfg.channel = setdiff(ftdata_posture.label,ft_channelselection({'eeg','O9','O10'},ftdata_posture));
        Z = ft_selectdata(cfg,ftdata_posture);
        X = regressOutByTrial(X,Y1);
        X = regressOutByTrial(X,Y2);
        X = regressOutByTrial(X,Y3);

        cfg = [];
        ftdata_posture = ft_appenddata(cfg,X,Z);

        %% get peak tremor from PCA
        [~, peakFrq,~,~,~,accData] = computeTremorProperties(ftdata_posture,[],0);
        cfg = [];
        ftdata_posture = ft_appenddata(cfg,ftdata_posture,accData);

        % Get rid of really big movements
        stdtr = cellfun(@(x) std(x')',accData.trial,'UniformOutput',0);
        sel = find(abs(zscore([stdtr{:}]))<3);
        cfg = [];
        cfg.trials = intersect(sel,trsel);
        ftdata_posture = ft_selectdata(cfg,ftdata_posture);

        %% Compute Spectra
        if strcmp(subTab.Hand{1},'L')
            list = {'Supp_Motor_Area_L*','Supp_Motor_Area_R*','Frontal_Sup_Medial_R*','Cerebellum_6_L*','Precuneus_R*'};
        else
            list = {'Supp_Motor_Area_R*','Supp_Motor_Area_L*','Frontal_Sup_Medial_L*','Cerebellum_6_R*','Precuneus_L*'};
        end

        cxSave = []; fxSave = []; cxMagSave = []; pxSave = [];
        for lroi = 1:5
            cfg = [];
            cfg.channel = list{lroi};
            cfg.method       = 'pca';
            cfg.numcomponent = 1;
            tmp = ft_componentanalysis(cfg,vc_data);
            tmp = rmfield(tmp,'topo');
            tmp = rmfield(tmp,'topolabel');
            tmp = rmfield(tmp,'unmixing');

            % Select on coherence level.
            cx = []; cmag = []; fx = [];
            triallist = intersect(sel,trsel);

            for trn = 1:numel(triallist)
                tr = triallist(trn);
                [px(:,trn),fpx] = pwelch(tmp.trial{tr},256,128,256,512);
                [cx(:,trn),fx] = mscohere(accData.trial{tr},tmp.trial{tr},256,128,512,512);
                fxind = find(fx>=4 & fx<=12);
                cmag(trn) = max(cx(fxind,trn));
            end
            pxSave(:,:,lroi) = px;
            cxSave(:,:,lroi) = cx;
            fxSave(:,lroi) = fx;
            fpxSave(:,lroi) = fpx;
            cxMagSave(:,lroi) = cmag;
        end
        tmp = {cxSave,fxSave,cxMagSave,pxSave,fpxSave};
        saveExpData(R,sub{1},'DICSImages','',[sub{1} 'PosturalHold_tremorDICSRegOut_SpecData_'],tmp,'Tremor');
        cxSaveSub{subI} = tmp{1};
        fxSaveSub{subI} = tmp{2};
        cxMagSaveSub{subI} = tmp{3};
        pxSaveSub{subI} = tmp{4};
        fpxSaveSub{subI} = tmp{5};

        %% Intersect the channels
        cfg = [];
        cfg.channel = sourceDICComFilt{1}.avg.label;
        cfg.channel = [cfg.channel; refch];
        ftdata_posture = ft_selectdata(cfg,ftdata_posture);

        % Take out the mean
        cfg = [];
        cfg.demean = 'yes';
        ftdata_posture = ft_preprocessing(cfg,ftdata_posture);

        % Apply the filter
        bnddef{1} = [peakFrq-1.5 peakFrq+1.5];
        bnd = 1;
        comfilt = sourceDICComFilt{bnd}.avg; %.filter;
        source_POST = runBeamFormer(ftdata_posture,headmodel,leadfield,'no',[],'DICSREF',bnddef{bnd},R.import.type,{refch});
        source_POST.pos = template_grid.pos;

        % Continue analysis
        cfg            = [];
        cfg.parameter = 'avg.coh';
        source_int  = ft_sourceinterpolate(cfg, source_POST, template_mri);

        % Z-normalize
        edf = 1/(size(source_POST.cumtapcnt,1)-2); % effective dof
        source_int.cohz = (atanh(source_int.coh)-edf)/sqrt(edf);

        % Flip LR
        subTab = retrieveSubjectInfoTab(R,sub);
        source_int.coh = flipVolume(table2array(subTab(:,'Hand')),source_int.coh,source_int.dim);
        source_int.cohz = flipVolume(table2array(subTab(:,'Hand')),source_int.cohz,source_int.dim);

        atlas = ft_convert_coordsys(atlas, source_int.coordsys);

        h = figure(subI); clf
        cfg              = [];
        cfg.method       = 'ortho';
        cfg.funparameter = 'coh';
        %             cfg.maskparameter = 'mask';
        %             cfg.funcolorlim   = [0 0.0175]; % wITH REGRESSION
        % cfg.funcolorlim   = [-0.1 0.1];
        cfg.location      = [-46 -20 64];
        cfg.funcolormap   = cmap;
        cfg.atlas = atlas;
        cfg.figure = h;
        ft_sourceplot(cfg,source_int);

        roi_center = [-6 0 70];
        radius = 1000; % mm
        global_peak = findLocalSourcePeak(source_int, roi_center, radius,'cohz');
        radius = 50; % mm
        local_peak = findLocalSourcePeak(source_int, roi_center, radius,'cohz');
        maxLabel = atlas_lookup_PublicFT(atlas,local_peak.pos,'coordsys','spm');
        locSave(:,subI) = {local_peak.pos local_peak.value global_peak.pos global_peak.value peakFrq maxLabel{1}};

        % % %% Save the Volume
        % % statvol = [];
        % % statvol.anatomy = reshape(source_int.cohz,source_int.dim);
        % % statvol.dim = source_int.dim;
        % % statvol.transform = source_int.transform;
        % % % get path
        % % path = [R.path.datapath '\' R.path.expname '\Group\PosturalHold\' sub{1} '_TremorDICS'];
        % % ft_write_mri([path '.nii'],statvol,'dataformat','nifti_spm')
        pathsave = saveExpData(R,sub{1},'DICSImages','',[sub{1} 'PosturalHold_tremorDICSRegOutVC_'],source_int,'Tremor');
        sourcesave{subI} = source_int;
    else
        sourcesave{subI} = loadExpData(R,sub{1},'DICSImages','',[sub{1} 'PosturalHold_tremorDICSRegOutVC_'],'Tremor');
    end
end

for subI = 1:numel(sourcesave)
    coh = sourcesave{subI}.coh;
    ZCoh = coh;%(coh-nanmean(coh))./nanstd(coh);
    sourcesave{subI}.coh = reshape(sourcesave{subI}.coh,sourcesave{subI}.dim);
    sourcesave{subI}.cohz = reshape(sourcesave{subI}.cohz,sourcesave{subI}.dim);
    sourcesave{subI}.ZCOH = reshape(ZCoh,sourcesave{subI}.dim);
    hsize = [10 10 10]; %[8 20]; % ~ 4 Hz x 500 ms gaussian
    sig = 3;
    sourcesave{subI}.ZCOHSm = convolve3DGaussian(sourcesave{subI}.ZCOH ,hsize,sig);
    sourcesave{subI}.dim = [91 109 91];
    sourcesave{subI}.freq = 5;
end

if EEGSel == 1
    % EEG List
    sourceCONT = sourcesave(1:11);
    sourceET = sourcesave([12:14 16:23]);
elseif EEGSel == 0
    %OPM List
    sourceCONT = sourcesave([1:4 8]);
    sourceET = sourcesave([5:7 9]);
end

%% Compute Overlap Figure
olap = []; subSig = [];
for subI = 1:numel(sourceET)
    olap(:,:,:,subI) = sourceET{subI}.ZCOHSm>prctile(sourceET{subI}.ZCOHSm(:),90);

    atlasLabel1 = find(strcmp(atlas.tissuelabel,'Cerebellum_Crus1_R'));
    atlasLabel2 = find(strcmp(atlas.tissuelabel,'Cerebellum_6_R'));
    x = [find(atlas.tissue == atlasLabel1); find(atlas.tissue == atlasLabel2)] ;
    xtmp = olap(:,:,:,subI);
    subSig(subI,1) = sum(xtmp(x)>prctile(sourceET{subI}.ZCOHSm(:),85))>(0.1*numel(x));

    atlasLabel1 = find(strcmp(atlas.tissuelabel,'Frontal_Sup_Medial_L'));
    x = find(atlas.tissue == atlasLabel1);
    xtmp = olap(:,:,:,subI);
    subSig(subI,2) = sum(xtmp(x)>prctile(sourceET{subI}.ZCOHSm(:),85))>(0.1*numel(x));

    atlasLabel1 = find(strcmp(atlas.tissuelabel,'Precuneus_L'));
    x = find(atlas.tissue == atlasLabel1);
    xtmp = olap(:,:,:,subI);
    subSig(subI,3) = sum(xtmp(x)>prctile(sourceET{subI}.ZCOHSm(:),85))>(0.1*numel(x));
end

olap = sum(olap,4);

overlapSource = sourcesave{1};
overlapSource.olap = olap;
overlapSource.threshmap = olap/max(olap(:));
cmapSet = brewermap(12,'RdBu');

h = figure(400+1); clf
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'olap';
cfg.funcolormap   = cmapSet;
% cfg.funcolorlim = [0 numel(sourceET)];
cfg.maskparameter    = 'threshmap';
cfg.maskstyle     = 'opacity';
cfg.atlas = atlas;
cfg.figure = h;
cfg.crosshair     = 'no';
ft_sourceplot(cfg,overlapSource);

%% Statistics below
% Plot T-statistic
stat = computeSourceTStatPaired(sourceET(find(subSig(:,1))),sourceCONT,'cohz'); % Select for coherence

hsize = [10 10 10]; %[8 20]; % ~ 4 Hz x 500 ms gaussian
sig = 3;
stat.stat = convolve3DGaussian(stat.stat,hsize,sig);

df = numel(sourceET)+numel(sourceCONT)-2;
stat.statmask = abs(stat.stat)>tinv(0.95,df);
stat.statmask(isnan(stat.statmask)) = 0;

h = figure(202); clf
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'stat';
cfg.maskparameter    = 'statmask';
cfg.maskstyle     = 'opacity';
% cfg.opacitylim    = [0 1];
% cfg.funcolorlim   = [-3 3];
cfg.funcolormap   = cmap;
% cfg.opacitymap    = 'vup';
% cfg.location = loclist{1};
cfg.atlas = atlas;
cfg.figure = h;
cfg.crosshair     = 'no';
cfg.axis = 'off';
ft_sourceplot(cfg,stat);

a = 1;
function signal = regressOutByTrial(signal,noise)
% Removes signal Y from X;
for tr = 1:numel(signal.trial)
    Sigtr = signal.trial{tr}; % Signal
    Noisetr = noise.trial{tr}; % nuisance
    for ch = 1:size(Sigtr,1)
        [~,res] = regressSignals(Noisetr',Sigtr(ch,:)',0);
        signal.trial{tr}(ch,:) = res';
    end
end
function zcoh = Znormcoh(coh)
edf = 1/(size(coh,2)-2); % effective dof
zcoh = (atanh(coh)-edf)/sqrt(edf);

