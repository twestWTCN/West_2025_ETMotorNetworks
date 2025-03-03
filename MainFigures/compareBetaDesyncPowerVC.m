function compareBetaDesyncPowerVC(R,sublist)
% Analyze beta desynchronization power in EEG/MEG channels during different task phases
% This function processes data to assess beta desynchronization across specific brain regions for different task parts, and visualizes results.

set(0,'defaultAxesFontSize',18); close all;

% List of channels to select for analysis
lrCh{1,1} = {'Supp_Motor_Area_L_*'}; % left hemisphere
lrCh{2,1} = {'Supp_Motor_Area_R_*'}; % right hemisphere

lrCh{1,2} = {'Precuneus_L_*'};
lrCh{2,2} = {'Precuneus_R_*'};

lrCh{1,3} = {'Frontal_Sup_Medial_L_*'};
lrCh{2,3} = {'Frontal_Sup_Medial_R_*'};

lrCh{1,4} = {'Cerebellum_6_L_*','Cerebellum_4_5_L_*'};
lrCh{2,4} = {'Cerebellum_6_R_*','Cerebellum_4_5_R_*'};

% Define time bounds for each phase
xlimz{2} = [-1.5 1.5];
xlimz{3} = [-1.5 1.5];
xlimz{4} = [-1.5 1.5];
xlimz{5} = [-1.5 3];

fresh = 1;
if fresh == 1
    subI = 0;
    for sub = sublist
        subI = subI + 1;
        for part = [2 3 4 5]
            % Setup data retrieval structure for analysis
            retdat.sub = sub;
            retdat.block = 3; % task block
            retdat.cond = 1:4; % condition indices
            retdat.part = part; % task phase
            retdat.ft = 2; % feature type
            retdat.ftflag = 1; % non-FT concatenation
            retdat.fileapend = [R.epoch.names{part} '_trans_pp_VC'];
            retdat.subfold = ['Condition' num2str(cond) '_epoched'];
            ftdata_cat = retrieveData(R, retdat);
            retdat.fileapend = [R.epoch.names{part} '_highRes_avtrl_VC_trans_freq'];
            freq_cat = retrieveData(R, retdat);

            % Temporarily remove 'pos' field and calculate descriptive statistics
            pos = freq_cat.pos;
            freq_cat = rmfield(freq_cat, 'pos');
            ftdata_cat = rmfield(ftdata_cat, 'pos');
            freq_cat = ft_freqdescriptives([], freq_cat);

            % Remove out-of-head channels
            cfg = [];
            cfg.channel = setdiff(freq_cat.label, ft_channelselection('oob_*', freq_cat.label));
            cfg.channel = setdiff(cfg.channel, ft_channelselection('*Pos', freq_cat.label));
            freq_cat = ft_selectdata(cfg, freq_cat);
            freq = freq_cat;

            % Normalize spectrogram by baseline
            if part == 2
                [freqBL, meanBase] = normalizeSpectrogramByBaseline(getPrePostEpochDetails(part), freq, 'relative');
            else
                [freqBL] = normalizeSpectrogramByBaseline([], freq, 'relative', meanBase);
            end

            % Flip sensors if subject is left-handed
            subTab = retrieveSubjectInfoTab(R, sub);
            if subTab.Hand{1} == 'L'
                freqBL = flipSensors(freqBL);
            end

            % Analyze selected channels
            for ch = 1:4
                for lr = 1:2
                    % Define time window for spectral power calculation
                    if part == 2 || part == 4
                        lat = getPrePostEpochDetails(part);
                    else
                        [~, lat] = getPrePostEpochDetails(part);
                    end

                    % Calculate signal-to-noise ratio (SNR) for each band
                    cfg = [];
                    cfg.latency = lat;
                    cfg.channel = lrCh{lr, ch};
                    ftdata_snr = ft_selectdata(cfg, ftdata_cat);

                    cfg = [];
                    cfg.method = 'mtmfft';
                    cfg.taper = 'dpss';
                    cfg.foilim = [2 48];
                    cfg.tapsmofrq = 1.5;
                    cfg.keeptrials = 'no';
                    freq_snr = ft_freqanalysis(cfg, ftdata_snr);
                    fx = squeeze(mean(freq_snr.powspctrm, 1));
                    hz = freq_snr.freq;

                    bndDef = [14 21; 21 30; 3 12];
                    for band = 1:3
                        betaInds = find(hz > bndDef(band, 1) & hz < bndDef(band, 2));
                        tailInds = find(hz > 30 & hz < 48);
                        snr(band, ch, lr, part, subI) = 10 .* log10(sum(fx(betaInds).^2) / sum(fx(tailInds).^2));
                    end

                    % Calculate desynchronization envelopes
                    cfg = [];
                    cfg.channel = lrCh{lr, ch};
                    freqLR = ft_selectdata(cfg, freqBL);

                    [pre, post] = getPrePostEpochDetails(part);
                    timeind = find(freqLR.time > post(1) & freqLR.time < post(2));

                    bndDef = [14 21; 21 30; 14 30];
                    for band = 1:3
                        bandind = find(freqLR.freq > bndDef(band, 1) & freqLR.freq < bndDef(band, 2));
                        tmp = mean(mean(freqLR.powspctrm(:, bandind, timeind), 2), 3);
                        dsTrace(:, ch, band, lr, part, subI) = squeeze(mean(mean(freqLR.powspctrm(:, bandind, :), 1), 2))';

                        if part == 3 || part == 5
                            [maxDesync(band, ch, lr, part, subI), ind] = max(tmp);
                        else
                            [maxDesync(band, ch, lr, part, subI), ind] = min(tmp);
                        end
                        meanDesync(band, ch, lr, part, subI) = median(tmp);
                    end
                end
            end
        end
    end
    % Save data based on data type
    if ~strcmp(R.import.type, 'OPM')
        save('tmp_VC_Deysncs'); % EEG
    else
        save('tmp_VC_deysncs_OPM'); % OPM
    end
else
    % Load precomputed data
    if ~strcmp(R.import.type, 'OPM')
        load('tmp_VC_Deysncs'); % EEG
    else
        load('tmp_VC_deysncs_OPM'); % OPM
    end
end

% Plot beta band desynchronization
figure(134); clf;
cmap = brewermap(4, 'Set1');
for band = 3
    for part = 2:5
        for lr = 1:2
            subplot(3, 4, sub2ind([4 3], part - 1, lr + 1));
            for ET = 1:2
                if ET == 1
                    if strcmp(R.import.type, 'OPM')
                        subInd{1} = [5 6 7 9];
                    else
                        subInd{1} = 1:11; % EEG
                    end
                elseif ET == 2
                    if strcmp(R.import.type, 'OPM')
                        subInd{2} = [2 3 4 8];
                    else
                        subInd{2} = [12:16 18:23];
                    end
                end
                XC = squeeze(dsTrace(:, ch, band, lr, part, subInd{ET}));
                subInd{ET} = subInd{ET}(max(abs(zscore(XC, [], 2))) < 3);
                XC = squeeze(dsTrace(:, ch, band, lr, part, subInd{ET}));
                XCSave{ET} = XC;

                % Plot mean and standard error
                XM = median(XC, 2);
                XSTD = std(XC, [], 2) ./ sqrt(size(XC, 2));
                [b, p(ET)] = boundedline(freqLR.time, XM, XSTD);
                b.Color = cmap(ET, :);
                b.LineWidth = 1.5;
                p(ET).FaceColor = cmap(ET, :);
                p(ET).FaceAlpha = 0.3;
                hold on;
                ylim([0.1 1.2]); xlim(xlimz{part});
                box off;
                axis square;
                grid on;
            end

            % Cluster-based permutation test
            dependent_samples = false;
            p_threshold = 0.1;
            num_permutations = 1000;
            two_sided = true;
            nC = 2;

            [clusters, p_values, t_sums, permutation_distribution] = permutest(XCSave{1}, XCSave{2}, dependent_samples, p_threshold, num_permutations, two_sided, nC);
            for i = 1:length(clusters)
                if p_values(i) < 0.05
                    line(freqLR.time(clusters{i}), repmat(1.1, size(clusters{i})), 'Color', cmap(2, :), 'LineWidth', 2);
                    text(freqLR.time(clusters{i}(1)), repmat(1.1, size(clusters{i}(1))), sprintf('%.3f', p_values(i)), 'Color', cmap(2, :), 'FontSize', 10);
                end
            end
            a = gca;
            a.Color = 'none';
        end
    end
    set(gcf, 'Position', [229 374 1024 512]);
end