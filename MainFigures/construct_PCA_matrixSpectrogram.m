function construct_PCA_matrixSpectrogram(R, sublist)
% Set some configurations
close all;
hrInd = 1;
hrtag = {'lowRes', 'highRes'};
partit = {'', 'Rest2Posture', 'Posture2Prep', 'Prep2Exec', 'Exec2Hold'};

% Recompute the main matrix of spectrograms?
fresh = 1;

% Initialize data containers
subI = 0;
freqMatCat = []; % Stores time-frequency data
trialInfoCat = []; % Stores trial information

% Loop through subjects and task parts
for sub = sublist
    subI = subI + 1;
    for part = [2 3 4 5]
        if fresh
            condpairs = {[1 3], [2 4], [1 2], [3 4]}; % Condition definitions - HUC, LUC, LRG, SML
            freq_cond = {};
            for cond = 1:4
                % Retrieve data for the specific condition and task part
                retdat.sub = sub;
                retdat.block = 3; % Task block
                retdat.cond = condpairs{cond}; % Condition
                retdat.part = part; % Movement phase
                retdat.ftflag = 1; % Non-FT concatenation
                retdat.subfold = ['Condition' num2str(cond) '_epoched'];
                
                % Specify data type (low or high resolution)
                switch hrtag{hrInd}
                    case 'lowRes'
                        retdat.fileapend = [R.epoch.names{part} '_lowRes_avtrl_VC_trans_freq'];
                    case 'highRes'
                        retdat.fileapend = [R.epoch.names{part} '_highRes_avtrl_VC_trans_freq'];
                end
                freq_cat = retrieveData(R, retdat);
                pos = freq_cat.pos;
                freq_cat = rmfield(freq_cat, 'pos');

                % Compute descriptive statistics per condition
                cfg = [];
                freq_cond{cond} = ft_freqdescriptives(cfg, freq_cat);
            end
            
            % Append data across conditions
            cfg = [];
            freq_cat = ft_appendfreq(cfg, freq_cond{:});

            % Identify patient or control and flip channels for left-handed subjects
            subTab = retrieveSubjectInfoTab(R, sub);
            if subTab.Hand{1} == 'L'
                freq_cat = flipSensors(freq_cat);
            end
            ET = cellfun(@(x) find(strcmp(x, {'CONT', 'ET'})), subTab.Condition) - 1;

            % Normalize by baseline
            freqBL = normalizeSpectrogramByBaseline(getPrePostEpochDetails(part), freq_cat, 'relative');

            % Select VC channels
            stringsToFind = {"Parietal_Sup_", "Precuneus_", "Supp_Motor_Area_", "Precentral_", "Postcentral_", "Frontal_Sup_", "Frontal_Mid_", "Cerebellum_6_", "Cerebellum_4_5_"};
            indexes = cell(1, numel(stringsToFind));
            for i = 1:numel(stringsToFind)
                searchString = stringsToFind{i};
                indexes{i} = find(contains(freqBL.label, searchString));
            end
            labList = freqBL.label(vertcat(indexes{:}));

            % Extract epochs of interest
            cfg = [];
            cfg.latency = [-1.5 3];
            cfg.frequency = [4 78];
            cfg.channel = labList;
            freqBL = ft_selectdata(cfg, freqBL);

            % Store positional data for channels
            posSel = pos(vertcat(indexes{:}), :);

            % Construct trial information
            trialInfoTmp = [1:4; repmat(subI, 4, 1)'; repmat(ET, 4, 1)'];

            % Preprocess the spectrogram data
            freqMatTmp = freqBL.powspctrm;
            freqMatTmp = log10(freqMatTmp);
            freqMatTmp(isnan(freqMatTmp)) = 0;
            freqMatTmp = (freqMatTmp - mean(freqMatTmp(:))) / std(freqMatTmp(:));

            % Save preprocessed data for further analysis
            saveExpData(R, sub{1}, 'PCAProcessedSpectrograms', '', [partit{part} '_PCAtrialinfo_' hrtag{hrInd}], trialInfoTmp, partit{part});
            saveExpData(R, sub{1}, 'PCAProcessedSpectrograms', '', [partit{part} '_PCAspectrograms_' hrtag{hrInd}], freqMatTmp, partit{part});
        else
            % Load precomputed spectrogram and trial information
            freqMatTmp = loadExpData(R, sub{1}, 'PCAProcessedSpectrograms', '', [partit{part} '_PCAspectrograms_' hrtag{hrInd}], partit{part});
            trialInfoTmp = loadExpData(R, sub{1}, 'PCAProcessedSpectrograms', '', [partit{part} '_PCAtrialinfo_' hrtag{hrInd}], partit{part});
        end
        
        % Compile data for PCA analysis
        freqMatCat(:,:,:,:,part-1,subI) = freqMatTmp; % cond x chan x freq x time x sub
        trialInfoCat(:,:,part-1,subI) = trialInfoTmp;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN ANALYSIS BEGINS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fresh
    save(['PosSEl_' hrtag{hrInd}], 'freqBL', 'posSel');
else
    load(['PosSEl_' hrtag{hrInd}]);
end
posSelScl = posSel * 10;

% Select data for group analysis and make a large matrix
if strcmp(R.import.type, 'EEG')
    freqMatCatX = freqMatCat(:,:,:,:,[2 3 4 5]-1, setdiff(1:23, [])); % cond x chan x freq x time x part x sub
    trialInfoCat = trialInfoCat(:,:,:, setdiff(1:23, []));
else
    freqMatCatX = freqMatCat(:,:,:,:,[2 3 4 5]-1, setdiff(1:9, [])); % cond x chan x freq x time x part x sub
    trialInfoCat = trialInfoCat(:,:,:, setdiff(1:9, []));
end
cmap = brewermap(4, 'Set1');
cmapSpectra = brewermap(256, '*RdYlBu');

% Plot group averages for each condition
figure(102); clf;
for part = 1:4
    for cond = 1:4
        subplot(4, 4, sub2ind([4 4], cond, part));
        imagesc(freqBL.time * 1000, freqBL.freq, (squeeze(mean(freqMatCatX(cond, :, :, :, part, :), [2 6]))));
        set(gca, 'YDir', 'normal');
        caxis([-0.6 0.6]);
        colormap(cmapSpectra);
        ylabel('Frequency (Hz)');
        xlabel('Time (ms)');
        axis square;
    end
end

%% Latent State Analysis Begins Here
XD = freqMatCatX;
XD(isinf(XD)) = 0;
[dim_cond, dim_chan, dim_freq, dim_time, dim_part, dim_sub] = size(XD);

% Reorder dimensions for PCA
XD = permute(XD, [3 2 4 5 1 6]); % Reorder: [freq chan time part cond sub]
N = dim_freq * dim_chan; % Parameters
SD = dim_time * dim_part; % Objects
REP = dim_cond * dim_sub; % Repeats

flattenedMatTr = reshape(XD, [N, SD, REP]); % [freq chan] x [time part] x [cond sub]
flattenedMat = mean(flattenedMatTr, 3); % Average across repeats
flattenedMat = flattenedMat'; % Transpose: [time part] x [freq chan]

%% PCA
NComp = 6; % Number of components
[score, mapping] = pca(flattenedMat, NComp);
coeff = mapping.M;
Rcoeff = rotatefactors(coeff);
mapping.name = 'PCA';

%% Scree Plot
latent = mapping.lambda;
latent = 100 .* (latent ./ sum(latent));

%% Save PCA Coefficients
coeffStruc = [];
coeffStruc.posSel = posSel;
coeffStruc.chLabel = freqBL.label;
coeffStruc.lambda = latent;
coeffStruc.score = score;
coeffStruc.coeff = coeff;
coeffStruc.Rotcoeff = Rcoeff;
coeffStruc.BP = coeff * coeff';
coeffStruc.dimList = '[freq chan] x [time part]';
coeffStruc.dataDim = [dim_cond, dim_chan, dim_freq, dim_time, dim_part, dim_sub];
coeffStruc.dataDimName = {'dim_cond', 'dim_channels', 'dim_freq', 'dim_time', 'dim_part', 'dim_sub'};
if strcmp(R.import.type, 'EEG')
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', coeffStruc, ['Coeff_' hrtag{hrInd}]);
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', freqMatCatX, ['freqMatCatX_' hrtag{hrInd}]);
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', trialInfoCat, ['trialInfoCat_' hrtag{hrInd}]);
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', freqBL, ['exampFreqBL_' hrtag{hrInd}]);
else
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', coeffStruc, ['Coeff_OPM_' hrtag{hrInd}]);
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', freqMatCatX, ['freqMatCatX_OPM_' hrtag{hrInd}]);
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', trialInfoCat, ['trialInfoCat_OPM_' hrtag{hrInd}]);
    saveExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', freqBL, ['exampFreqBL_OPM_' hrtag{hrInd}]);
end
