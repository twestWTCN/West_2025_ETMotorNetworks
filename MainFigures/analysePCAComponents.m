function analysePCAComponents(R)
% Analyze PCA Components from spectrogram data
% This function analyzes the principal components derived from spectrogram data, plotting both spatial and spectral characteristics.

cmap = brewermap(6, 'Set1');
cmapSpectra = brewermap(256, '*RdYlBu');
lspec = {'-', '--', ':', '-.'};
partit = {'', 'Rest2Posture', 'Posture2Prep', 'Prep2Exec', 'Exec2Hold'};
hrtag = {'lowRes', 'highRes'};

% Load in the precomputed coefficients
hrInd = 1; % High-resolution version
mesh2 = load(fullfile(R.path.ftpath, 'template/anatomy/surface_pial_both'), 'mesh'); % Load mesh

if strcmp(R.import.type, 'EEG')
    coeffStruc = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['Coeff_' hrtag{hrInd}]);
    freqMatCatX = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['freqMatCatX_' hrtag{hrInd}]);
    trialInfoCat = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['trialInfoCat_' hrtag{hrInd}]);
    freqBL = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['exampFreqBL_' hrtag{hrInd}]);
else
    coeffStruc = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['Coeff_' hrtag{hrInd}]);
    freqMatCatX = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['freqMatCatX_OPM_' hrtag{hrInd}]);
    trialInfoCat = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['trialInfoCat_OPM_' hrtag{hrInd}]);
    freqBL = loadExpData(R, 'GROUP', 'PCA', [], 'Coefficients&Dims', ['exampFreqBL_OPM_' hrtag{hrInd}]);
    coeffStruc.dataDim(6) = 9; % Adjust subject number for OPM
end

xlimz{2} = [-1.5 1.5];
xlimz{3} = [-1.5 1.5];
xlimz{4} = [-1.5 1.5];
xlimz{5} = [-1.5 3];

%% Unpack Data Matrix
XD = freqMatCatX;
XD(isinf(XD)) = 0;
[dim_cond, dim_chan, dim_freq, dim_time, dim_part, dim_sub] = size(XD);

% Reorder dimensions for analysis
XD = permute(XD, [3 2 4 5 1 6]); % Reorder: [freq chan time part cond sub]
N = dim_freq * dim_chan; % Parameters
SD = dim_time * dim_part; % Objects
REP = dim_cond * dim_sub; % Repeats

flattenedMatTr = reshape(XD, [N, SD, REP]); % [freq chan] x [time part] x [cond sub]
flattenedMat = mean(flattenedMatTr, 3); % Average across repeats
flattenedMat = flattenedMat'; % Transpose: [time part] x [freq chan]

%% Unpack the Coefficients and Scores
NComp = 6;
coeff = coeffStruc.Rotcoeff;
score = coeffStruc.score;
posSel = coeffStruc.posSel * 10;

coeff = coeff(:, 1:NComp);
score = score(:, 1:NComp);

reshapeCoeff = reshape(coeff, [dim_freq dim_chan NComp]); % [hz x chan x comp]
reshapeScore = reshape(score, [dim_time dim_part NComp]); % [time x part x comp]
latentTime = reshapeScore;

%% Plot Each State's Parameters (Frequency and Space)
% Get peak frequency of each channel and each state
coeffSpectra_M = squeeze(mean(reshapeCoeff, [2])); % Mean spectra across channels
coeffSpectra_S = squeeze(std(reshapeCoeff, [], [2])); % Std spectra across channels

% Distribution of peak frequencies within channels
for comp = 1:4
    for ch = 1:dim_chan
        [~, pki] = max(abs(squeeze(reshapeCoeff(:, ch, comp))));
        peakFrqCh(ch, comp) = freqBL.freq(pki);
        peakFrqLoc(ch, comp) = pki;
    end
end

peakFrqLocTrunc = peakFrqLoc;
peakFrqLocTrunc(peakFrqLoc > 48) = 48;
cmapSpectraFlat = brewermap(48, '*RdYlBu');

% Get Channel Weights
betaInd{1} = freqBL.freq >= 14 & freqBL.freq <= 21;
coeffChannelweight{1} = squeeze(sum(reshapeCoeff(betaInd{1}, :, 1), 1));

betaInd{2} = freqBL.freq >= 30 & freqBL.freq <= 60;
coeffChannelweight{2} = squeeze(sum(reshapeCoeff(betaInd{2}, :, 2), 1));

betaInd{3} = freqBL.freq >= 8 & freqBL.freq <= 15;
coeffChannelweight{3} = squeeze(sum(reshapeCoeff(betaInd{3}, :, 3), 1));

betaInd{4} = freqBL.freq >= 18 & freqBL.freq <= 34;
coeffChannelweight{4} = squeeze(sum(reshapeCoeff(betaInd{4}, :, 4), 1));

%% Plotting Components
figure(231); clf;
for comp = 1:4
    % Plot the spectra
    subplot(6, 6, sub2ind([6 6], comp + 1, 1));
    [bspec, pspec] = boundedline(freqBL.freq, coeffSpectra_M(:, comp), coeffSpectra_S(:, comp));
    bspec.Color = cmap(comp, :); bspec.LineWidth = 1.5;
    pspec.FaceColor = cmap(comp, :); pspec.FaceAlpha = 0.4;
    hold on;
    axis square; grid off; box off;
    ylim([-0.025 0.035]);

    % Plot the histograms of peak frequencies
    yyaxis right;
    h(comp) = histogram(peakFrqCh(:, comp), [0:2:68], 'Normalization', 'probability');
    h(comp).FaceColor = cmap(comp, :);
    xlim([0 48]);
    a = gca;
    a.XTickLabel = [];
    a.YColor = 'k';

    % Plot the spatial distribution
    subplot(6, 6, sub2ind([6 6], comp + 1, 2));
    brain = mesh2.mesh;
    hs = ft_plot_mesh(brain);
    camlight;
    hs.FaceAlpha = 0.05;
    hold on;

    % Adjust and threshold sizes
    Sz = coeffChannelweight{comp};
    if comp ~= 2
        SzInd = Sz > prctile(Sz, 85);
    else
        SzInd = Sz < prctile(Sz, 15);
    end
    Sz = Sz(SzInd);

    scatter3(posSel(SzInd, 1), posSel(SzInd, 2), posSel(SzInd, 3), 100 .* (abs(Sz) ./ max(abs(Sz))), cmap(comp, :), 'filled', 'MarkerEdgeColor', 'none');
    view(90, 0);
end

% Plot empirical and reconstructed data for each part
for part = 1:4
    subplot(6, 6, sub2ind([6 6], 1, part + 2));
    imagesc(freqBL.time, freqBL.freq, (squeeze(mean(freqMatCatX(:, :, :, :, part, :), [1 2 6]))));
    set(gca, 'YDir', 'normal');
    colormap(cmapSpectra);
    axis square; box off;
    xlim(xlimz{part + 1});

    for comp = 1:4
        subplot(6, 6, sub2ind([6 6], comp + 1, part + 2));
        BP = coeff(:, comp) * coeff(:, comp)';
        X = flattenedMat;
        X = X * BP;
        X = reshape(X, [dim_time dim_part dim_freq dim_chan]);
        imagesc(freqBL.time, freqBL.freq, squeeze(mean(X(:, part, :, :), 4))');
        set(gca, 'YDir', 'normal'); box off;
        xlim(xlimz{part + 1});
        colormap(cmapSpectra);
        axis square;
    end

    subplot(6, 6, sub2ind([6 6], 6, part + 2));
    BP = coeff * coeff';
    X = flattenedMat;
    X = X * BP;
    X = reshape(X, [dim_time dim_part dim_freq dim_chan]);
    imagesc(freqBL.time, freqBL.freq, squeeze(mean(X(:, part, :, :), 4))');
    set(gca, 'YDir', 'normal'); box off;
    xlim(xlimz{part + 1});
    colormap(cmapSpectra);
    axis square;
end
set(gcf, 'Position', [400 100 1400 800]);

%% Render Spatial Map again separately
figure(821)
for comp = 1:4
    i = 1;
    subplot(1, 4, comp)
    brain = mesh2.mesh;
    hs = ft_plot_mesh(brain);
    camlight;
    hs.FaceAlpha = 0.05;
    hold on;

    % Adjust and threshold sizes
    Sz = coeffChannelweight{comp};
    if comp ~= 2
        SzInd = Sz > prctile(Sz, 85);
    else
        SzInd = Sz < prctile(Sz, 15);
    end
    Sz = Sz(SzInd);

    scatter3(posSel(SzInd, 1), posSel(SzInd, 2), posSel(SzInd, 3), 100 .* (abs(Sz) ./ max(abs(Sz))), cmap(comp, :), 'filled', 'MarkerEdgeColor', 'none');
    axis equal;
    if i == 1; view(90, 0); elseif i == 2; view(0, 90); end
end
set(gcf, 'Position', [400 100 1800 900])

%% Project Single Trial Latents
trialLatents = [];
for comp = 1:NComp
    for i = 1:size(flattenedMatTr, 3)
        XTmp = flattenedMatTr(:, :, i); % [freq chan] x [time part] x [cond sub]
        XTmp = XTmp - mean(XTmp, 1);
        XTmp = XTmp' * coeff(:, comp);
        XTmp = XTmp + mean(flattenedMatTr(:, :, i), 1)';
        XTmp = reshape(XTmp, [dim_time dim_part]);
        trialLatents(:, :, i, comp) = XTmp; % time x part x [dim_cond dim_sub] x comp
    end
end

%% Plotting Latents Per Condition
XStore1 = reshape(trialLatents, [dim_time dim_part dim_cond dim_sub NComp]);
figure(121); clf;
for comp = 1:4
    for part = 1:4
        subplot(4, 4, sub2ind([4 4], comp, part))
        for cond = 1:4
            [bl(cond), bp] = boundedline(freqBL.time * 1000, mean(XStore1(:, part, cond, :, comp), 4), std(XStore1(:, part, cond, :, comp), [], 4) ./ sqrt(size(XStore1, 4)));
            bl(cond).Color = cmap(cond, :);
            bl(cond).LineWidth = 1.5;
            bp.FaceColor = cmap(cond, :);
            bp.FaceAlpha = 0.1;
            hold on;
        end
        ylabel(sprintf('Comp %.0f', comp))
        title(R.epoch.names{part + 1})
    end
end

%% Plotting Latents Per Patient vs Control
ETSub = squeeze(trialInfoCat(3, 1, 1, :)) + 1;
p_valuesS = []; pbank = [];
figure(122); clf;
for comp = 1:4
    partI = 0;
    for part = 1:4
        partI = partI + 1;
        subplot(6, 4, sub2ind([4 6], partI, comp + 1))
        for ET = 1:2
            tmp = squeeze(XStore1(:, part, :, ETSub == ET, comp));
            XCSave{ET} = reshape(tmp, size(tmp, 1), prod(size(tmp, [2 3])));

            % Weights for gamma reversed
            if comp == 2
                invt = -1;
            else
                invt = 1;
            end

            [bl(ET), bp] = boundedline(freqBL.time, invt * mean(XStore1(:, part, :, ETSub == ET, comp), [3 4]), invt * std(XStore1(:, part, :, ETSub == ET, comp), [], [3 4]) ./ sqrt(sum(ETSub == ET)));
            bl(ET).Color = cmap(ET, :);
            bl(ET).LineWidth = 1.5;
            bp.FaceColor = cmap(ET, :);
            bp.FaceAlpha = 0.3;
            hold on;
        end

        % Perform permutation test
        dependent_samples = false; % Assuming the samples are dependent
        p_threshold = 0.1; % Significance level
        num_permutations = 5000; % Number of permutations
        two_sided = true; % Two-sided test
        nC = 2;
        [clusters, p_values, t_sums, permutation_distribution] = permutest(XCSave{1}, XCSave{2}, dependent_samples, p_threshold, num_permutations, two_sided, nC);
        for i = 1:length(clusters)
            if p_values(i) < 0.05
                line(freqBL.time(clusters{i}), repmat(20, size(clusters{i})), 'Color', 'k', 'LineWidth', 4);
            end
        end
        pbank = [pbank p_values];
        axis([-1500 3000 -50 50])
        xlim(xlimz{part + 1})
        axis on; axis square;
    end
end
set(gcf, 'Position', [672 230 953 656])

%% Scatter Plot of Components
[~, zeroind] = min((freqBL.time - 0) .^ 2);
sc = [];
partI = 0;
figure(123); clf;
for part = 1:4
    partI = partI + 1;
    subplot(6, 4, sub2ind([4 6], partI, 6))
    for ET = 1:2
        a = squeeze(mean(XStore1(:, part, :, ETSub == ET, 1), [3 4]));
        b = squeeze(mean(XStore1(:, part, :, ETSub == ET, 4), [3 4]));
        scatter(mean(a(1, :), 2), mean(b(1, :), 2), 100, cmap(ET, :), 'o', 'filled')
        hold on;
        a = smooth(mean(a, 2), 5); b = smooth(mean(b, 2), 5);
        sc(ET) = scatter(a(zeroind, :), b(zeroind, :), 100, cmap(ET, :), 'x', 'LineWidth', 2);
        scatter(a(end, :), b(end, :), 50, cmap(ET, :), '>', 'filled', 'LineWidth', 2);
        plot(a, b, 'Color', cmap(ET, :), 'LineWidth', 2);

        axis off; axis square;
    end
    if part == 4
        legend(sc, {'Cont', 'ET'}, 'NumColumns', 2, 'Box', 'off', 'Location', 'NorthEastOutside');
    end
end
set(gcf, 'Position', [700 100 900 800])

% Plot the traces
for part = 1:4
    subplot(6, 4, sub2ind([4 6], part, 1))
    plot(acctrace{part + 1}(1, :), acctrace{part + 1}(2, :), 'k')
    ylim([-0.3 0.3])
    xlim(xlimz{part + 1})
    box off;
    axis square;
end
