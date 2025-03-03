function compareTaskSpectrograms_VC(R,sublist)
% Compare spectrograms across different experimental task phases
% This function retrieves spectrogram data for specified subjects and experimental conditions,
% performs preprocessing including channel selection and baseline normalization, and then generates
% spectrogram visualizations with overlaid statistical analysis. It is intended for use with EEG or MEG data.

set(0,'defaultAxesFontSize',18); close all
partit = {'','Rest2Posture','Posture2Prep','Prep2Exec','Exec2Hold'};

fresh = 1;
subI = 0;
for sub = sublist
    subI = subI + 1;
    for part = 2:5
        if fresh == 1
            % Setup retrieval structure for task data
            retdat.sub = sub;
            retdat.block = 3; % Task block
            retdat.cond = 1:4; % Condition indices
            retdat.part = part; % Movement phases
            retdat.ft = 2; % Feature type
            retdat.ftflag = 1; % Feature flag for non-FT concatenation
            retdat.fileapend = [R.epoch.names{part} '_trans_pp_arV3'];
            retdat.subfold = ['Condition' num2str(cond) '_epoched'];
            ftdata_cat = retrieveData(R, retdat);
            retdat.fileapend = [R.epoch.names{part} '_highRes_avtrl_VC_trans_freq'];
            freq_cat = retrieveData(R, retdat);

            % Temporarily remove problematic field
            postmp = freq_cat.pos;
            freq_cat = rmfield(freq_cat, 'pos');
            freq_cat = ft_freqdescriptives([], freq_cat);

            % Get channel data for relevant signals
            cfg = [];
            pattern = '.*Finger1Pos';
            chInd = find(~cellfun('isempty', regexp(ftdata_cat.label, pattern)));
            cfg.channel = ftdata_cat.label(chInd);
            ftdata_catA = ft_selectdata(cfg, ftdata_cat);

            cfg = [];
            cfg.channel = {'Coder'};
            ftdata_catB = ft_selectdata(cfg, ftdata_cat);

            % Standardize trials for each subject
            for trn = 1:numel(ftdata_catA.trial)
                ftdata_catA.trial{trn} = ft_preproc_standardize(ftdata_catA.trial{trn});
                ftdata_catB.trial{trn} = ftdata_catB.trial{trn} * 20;
            end
            freq = freq_cat;

            % Perform baseline normalization
            cfg = [];
            cfg.baseline = getPrePostEpochDetails(part);
            cfg.baselinetype = 'zscore';
            freqBL = ft_freqbaseline(cfg, freq);

            % Flip channels if subject is left-handed
            subTab = retrieveSubjectInfoTab(R, sub);
            if subTab.Hand{1} == 'L'
                freq = flipSensors(freqBL);
            end

            % Select Channels of interest
            cfg = [];
            cfg.channel = {'Precuneus_L_*', 'Precuneus_R_*', 'Precentral_L_*', 'Precentral_R_*',
                           'Supp_Motor_Area_L_*', 'Supp_Motor_Area_R_*', 'Frontal_Mid_L*', 'Frontal_Sup_L_*',
                           'Frontal_Mid_R*', 'Frontal_Sup_R_*', 'Cerebellum_6_R_*', 'Cerebellum_6_L_*',
                           'Cerebellum_4_5_L_*', 'Cerebellum_4_5_R_*'};
            freqBL = ft_selectdata(cfg, freqBL);

            % Store results for each subject and part
            codmotStore{subI, part} = [mean(vertcat(ftdata_catA.trial{:})); mean(vertcat(ftdata_catB.trial{:}))];
            codmotStoreTVEC{subI, part} = mean(vertcat(ftdata_catA.time{:}));
            freqStore{subI, part} = freqBL;

            % Create spectrogram structure and save
            spect.freqBL = freqBL;
            spect.codmotStore = [mean(vertcat(ftdata_catA.trial{:})); mean(vertcat(ftdata_catB.trial{:}))];
            spect.codmotStoreTVEC = mean(vertcat(ftdata_catA.time{:}));
            saveExpData(R, sub{1}, 'VCSpectrogramsV3', '', [partit{part} '_spectrograms'], spect, partit{part});
        else
            % Load precomputed spectrogram data
            spect = loadExpData(R, sub{1}, 'VCSpectrogramsV3', '', [partit{part} '_spectrograms'], partit{part});
            codmotStore{subI, part} = spect.codmotStore;
            codmotStoreTVEC{subI, part} = spect.codmotStoreTVEC;
            freqStore{subI, part} = spect.freqBL;
        end
    end
end

% Set subjects selection for different groups
subSelOPMCont = {[], setdiff(1:5, []), setdiff(1:5, []), setdiff(1:5, []), setdiff(1:5, [])};
subSelOPMPat = {[], setdiff(1:4, []), setdiff(1:4, []), setdiff(1:4, []), setdiff(1:4, [])};

subSelEEGCont = {[], setdiff(1:11, []), setdiff(1:11, []), setdiff(1:11, []), setdiff(1:11, [])};
subSelEEGPat = {[], setdiff(1:12, []), setdiff(1:12, []), setdiff(1:12, []), setdiff(1:12, [])};

subSel = subSelOPMCont;

% Grand Average and Plotting
cmap = brewermap(256, '*RdYlBu');

figure(101); clf
chnames = {'Supp_Motor_Area_L_*', 'Precuneus_L_*', {'Cerebellum_6_R_*', 'Cerebellum_4_5_R_*'}};
for chtype = 1:3
    for part = 2:5
        XTmp = [];
        for i = 1:numel(subSel{part})
            XD = freqStore{subSel{part}(i), part};
            cfg = [];
            cfg.channel = chnames{chtype};
            XD = ft_selectdata(cfg, XD);
            d = squeeze(mean(XD.powspctrm, 1));
            XTmp(:, i) = d(:);
        end
        % T-statistics calculation and plotting
        XBar = mean(XTmp, 2);
        SE = std(XTmp, [], 2) ./ sqrt(numel(subSel{part}));
        T = reshape(XBar ./ SE, size(d));
        
        hsize = [64 120];
        sig = 1;
        TConv = convolve2DGaussianWithNans(T, hsize, sig);

        cfg = [];
        cfg.channel = chnames{chtype};
        freqGA = ft_freqgrandaverage(cfg, freqStore{subSel{part}, part});

        f = subplot(3, 4, sub2ind([4 3], part-1, chtype));
        imagesc(freqGA.time * 1000, freqGA.freq, TConv);
        set(gca, 'YDir', 'normal');
        hold on

        threshmask = abs(TConv) > 1.96;
        [C, h] = contour(freqGA.time * 1000, freqGA.freq, threshmask, [1 1], 'LineColor', 'k', 'LineWidth', 2);
        
        caxis([-5 5]);
        if part == 5
            C = colorbar('Location', 'EastOutside');
        end
        xlim([-1000 1500]); ylim([4 68]);
        xlabel('Time from onset (ms)');
        if part == 2
            ylabel('Frequency (Hz)');
        end
        axis square;
        colormap(cmap);

        % Plot kinematics
        comot = vertcat(codmotStore{subSel{part}, part});
        comot = comot(1:2:end, :);
        comot = abs(comot - comot(:, 1));
        gcmp = nanmean(comot);
        yyaxis right; axis square;
        plot(codmotStoreTVEC{1, part} * 1000, gcmp, 'b--', 'LineWidth', 2);
        ylim([0 1.5]);

        a = gca;
        a.YTickLabel = '';
    end
    yyaxis left;
    C.Position = [0.9182 0.3063 0.0130 0.4203];
    C.Label.String = 'T-score';
    set(gcf, 'Position', [211 43 1570 935]);
end
