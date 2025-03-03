function projectPCA2Trials(R, sublist)
%% Function projectPCA2Trials
% This function performs a PCA projection on frequency spectrogram data across different experimental conditions.
% It retrieves, normalizes, and processes trial data for each subject and condition before applying PCA.
% The script also extracts kinematic and tremor-related information from trial metadata.
% 
% Arguments:
% R - Experiment configuration and data information
% sublist - List of subjects to include in the analysis
% Define resolution tag (either low or high resolution)
hrtag = {'lowRes', 'highRes'};
hrInd = 1; % Index for selecting resolution version

% Define colormap for spectral plots
cmapSpectra = brewermap(256, '*RdYlBu');
cmapset = brewermap(12, 'Set1');
cmap = brewermap(256, '*RdYlBu');

% Close all open figures and set fresh flag to load new data
close all;
fresh = 1;

% Load PCA coefficients
coeffStruc = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['Coeff_' hrtag{hrInd}]);
coeff = coeffStruc.Rotcoeff;

% Loop through specific parts of the trial data
for part = 3:5
    if fresh
        subI = 0;
        subX = []; tremDetX = []; kinDetX = []; latentX = [];

        % Iterate over each subject in the subject list
        for sub = sublist
            subI = subI + 1;
            freq_cond = cell(1, 4); % Pre-allocate for four conditions

            % Iterate over each condition (1-4)
            for cond = 1:4
                % Set data retrieval parameters for each condition
                retdat.sub = sub;
                retdat.block = 3; % Task block
                retdat.cond = cond; % Experimental condition
                retdat.part = part; % Movement part
                retdat.ftflag = 1; % Non-FT concatenation
                retdat.subfold = ['Condition' num2str(cond) '_epoched'];

                % Append appropriate filename based on resolution
                switch hrtag{hrInd}
                    case 'lowRes'
                        retdat.fileapend = [R.epoch.names{part} '_lowRes_keeptrl_VC_trans_freq'];
                    case 'highRes'
                        retdat.fileapend = [R.epoch.names{part} '_highRes_keeptrl_VC_trans_freq'];
                end
                freq_cond{cond} = retrieveData(R, retdat);
                pos = freq_cond{cond}.pos;
                freq_cond{cond} = rmfield(freq_cond{cond}, 'pos');

                % Add trial info from additional data retrieval
                retdat.fileapend = [R.epoch.names{part} '_trans_pp_arV3'];
                dum_data = retrieveData(R, retdat);
                freq_cond{cond}.trialinfo = dum_data.trialinfo;
                clear dum_data;
            end

            % Concatenate conditions together
            cfg = [];
            freq_cat = ft_appendfreq(cfg, freq_cond{:});

            % Retrieve subject-specific information and handle handedness
            subTab = retrieveSubjectInfoTab(R, sub);
            if subTab.Hand{1} == 'L'
                freq_cat = flipSensors(freq_cat); % Flip channels for left-handed subjects
            end

            % Normalize spectrogram by baseline
            freqBL = normalizeSpectrogramByBaseline(getPrePostEpochDetails(part), freq_cat, 'relative', [], 1);

            % Select channels of interest (VC Channels)
            stringsToFind = {"Parietal_Sup_", "Precuneus_", "Supp_Motor_Area_", "Precentral_", "Postcentral_", "Frontal_Sup_", "Frontal_Mid_", "Cerebellum_6_", "Cerebellum_4_5_"};
            indexes = cell(1, numel(stringsToFind));
            for i = 1:numel(stringsToFind)
                searchString = stringsToFind{i};
                indexes{i} = find(contains(freqBL.label, searchString));
            end
            labList = freqBL.label(vertcat(indexes{:}));

            % Select data within desired latency and frequency range
            cfg = [];
            cfg.latency = [-1.5 3];
            cfg.frequency = [4 78];
            cfg.channel = labList;
            freqBL = ft_selectdata(cfg, freqBL);

            % Get positional data and assign trial info
            posSel = pos(vertcat(indexes{:}), :);
            labelList{subI} = labList;
            trialInfoTmp = freqBL.trialinfo;
            trialInfoTmp = [trialInfoTmp repmat(part, size(trialInfoTmp, 1), 1)];

            % Preprocess time-frequency matrix (TF)
            freqMatTmp = freqBL.powspctrm;
            freqMatTmp = log10(freqMatTmp);
            freqMatTmp(isnan(freqMatTmp)) = 0;

            % Smooth the individual TF trials using Gaussian filter
            [dim_rep, dim_chan, dim_freq, dim_time] = size(freqMatTmp);
            freqMatTmpSm = nan(size(freqMatTmp));
            for rep = 1:dim_rep
                for chan = 1:dim_chan
                    X = squeeze(freqMatTmp(rep, chan, :, :));
                    freqMatTmpSm(rep, chan, :, :) = convolve2DGaussianWithNans(X, [8 8], 1);
                end
            end
            % Normalize in the time domain
            freqMatTmpSm = (freqMatTmpSm - mean(freqMatTmpSm, 3)) ./ std(freqMatTmpSm, [], 3);

            % Reshape for PCA projection
            XD = permute(freqMatTmpSm, [3 2 4 1]); % Reorder: [freq chan time rep]
            N = dim_freq * dim_chan;
            SD = dim_time;
            REP = dim_rep;
            flattenedMatTr = reshape(XD, [N, SD, REP]);

            % Project the data using PCA coefficients
            backProjData = [];
            latent = [];
            for comp = 1:4
                for i = 1:size(flattenedMatTr, 3)
                    XTmp = flattenedMatTr(:, :, i);
                    XTmp = XTmp - mean(XTmp, 1);
                    XTmp = XTmp' * coeff(:, comp);
                    latent(:, i, comp) = XTmp + mean(flattenedMatTr(:, :, i), 1)';
                end
            end

            % Average across channels and reshape
            backProjData = squeeze(mean(backProjData, 1));
            backProjData = permute(backProjData, [3 4 1 2]); % rep x comp x freq x chan

            % Concatenate data for each subject
            subX = cat(1, subX, backProjData);
            latentX = cat(2, latentX, latent);

            % Extract tremor and kinematic details
            tremDet = trialInfoTmp(:, [7 15:18]);
            tremDetX = cat(1, tremDetX, tremDet);
            kinDet = trialInfoTmp(:, [7 18:end]);
            kinDetX = cat(1, kinDetX, kinDet);
        end

        % Save processed data
        save(['E:\TimWest\GITHUB\ReachingToUnderstandET_V2\tmpdata\PCA_DECODERDATA_' R.epoch.names{part} '_LowRes.mat'], "tremDetX", "kinDetX", "subX", "freqBL", "latentX", "parlabs")
    else
        % Load pre-processed data
        load(['E:\TimWest\GITHUB\ReachingToUnderstandET_V2\tmpdata\PCA_DECODERDATA_' R.epoch.names{part} '_LowRes.mat'], "tremDetX", "kinDetX", "subX", "freqBL", "latentX", "parlabs")
    end
end
end