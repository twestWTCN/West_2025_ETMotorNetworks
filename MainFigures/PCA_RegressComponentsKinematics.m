function PCA_RegressComponentsKinematics(R)
%% Function PCA_RegressComponentsKinematics
% This function performs regression analysis on PCA components against kinematic and tremor-related features.
% It loads preprocessed data, computes regressors, and performs statistical analysis and visualization to assess relationships.
%
% Arguments:
% R - Experiment configuration and data information

% Set colormap and figure parameters
cmapset = brewermap(14, 'Paired');
cmapsetE{1} = cmapset(1:2:end, :);
cmapsetE{2} = cmapset(2:2:end, :);
cmap = brewermap(256, '*RdYlBu');

% Define parameter labels
parTit = {'SubI', 'tremVar', 'mean_velocity', 'max_velocity', 'pathlength', 'mov_dur', 'smooth_index', 'trialDeviance', 'holdStability', 'reachDur', 'reachRT', 'holdRT', 'postremPow', 'partN'};
parTitTrem = {'SubI', 'tremPeakPow', 'tremFreq', 'tremPow', 'tremVar'};
regFeatPars = {'tremPeakPow', 'tremFreq', 'tremPow', 'tremVar', 'mean_velocity', 'pathlength', 'smooth_index', 'holdStability', 'reachDur', 'reachRT', 'holdRT'};

close all;
fresh = 0;
partI = 0;
savea = [];

% Iterate through parts of the trial data (e.g., different conditions)
for part = [2 3 4 5]
    partI = partI + 1;

    % Load preprocessed data
    load(['E:\TimWest\GITHUB\ReachingToUnderstandET_V2\tmpdata\PCA_DECODERDATA_' R.epoch.names{part} '_LowRes.mat'], "tremDetX", "latentX", "subX", "kinDetX", "freqBL", "parlabs");
    load('E:\TimWest\GITHUB\ReachingToUnderstandET_V2\tmpdata\MovementVecs.mat');

    % Define axis limits for different parts
    xlimz{2} = [-1.500 1.500];
    xlimz{3} = [-1.500 1.500];
    xlimz{4} = [-1.500 1.500];
    xlimz{5} = [-1.500 3.000];

    % Set up regression features based on the part of the data
    subInds = kinDetX(:, 1);
    ETG = ones(size(subInds));
    ETG(subInds > 11) = 2;

    if part == 2 || part == 3
        regressorFeat = [tremDetX(:, 2:end), kinDetX(:, [3 6 7 9 10])];
    else
        regressorFeat = [tremDetX(:, 2:end), kinDetX(:, [3 6 7 9 10 11 12])];
    end

    % Remove outliers based on z-score
    subCodes = unique(subInds);
    for i = 1:numel(subCodes)
        rowIdx = subInds == subCodes(i);
        XP = regressorFeat(rowIdx, :);
        outThresh = nanzscore(XP, 1);
        XP(abs(outThresh) > 3) = nan;
        regressorFeat(rowIdx, :) = XP;
    end

    % Define feature titles for plotting
    [~, zeroind] = min((freqBL.time - 0) .^ 2);
    featTit = [parTitTrem(2:end), parTit([3 5 7 9 10 11 12])];

    %% Run regression and visualization for each feature
    featI = 0;
    for featInd = [1 5 8] % Specific indices of regressor features
        featI = featI + 1;

        % Define subjects for analysis
        if featInd == 1
            sublist = [12:14, 16:23]; % Only compute tremor stats for ET subjects
        else
            sublist = 1:23;
        end

        figure(200);
        offset = -8 * (featI - 1);

        for comp = 1:4
            for condcomp = 1:2
                % Split data for each subject
                tmp = [];
                for i = unique(sublist)
                    L = regressorFeat(subInds == i, featInd);
                    tmpX = squeeze(latentX(:, subInds == i, comp));
                    if condcomp == 1
                        tmpX = tmpX(:, L <= prctile(L, 25));
                    else
                        tmpX = tmpX(:, L >= prctile(L, 75));
                    end
                    XC{condcomp, i} = tmpX;
                    tmp = [tmp, nanmean(tmpX, 2)];
                end

                % Compute mean and standard error
                X = 2 * squeeze(mean(tmp, 2)); % Scaling for visualization
                XS = 2 * squeeze(std(tmp, [], 2) ./ sqrt(size(tmp, 2))); % Standard error

                % Reverse weights for gamma component
                if comp == 2
                    X = -X;
                    XS = -XS;
                end

                % Plot bounded line
                XComp{condcomp, comp} = tmp;
                subplot(5, 4, sub2ind([4 5], partI, comp));
                [bl(condcomp), bp] = boundedline(freqBL.time, X + offset, XS);
                bl(condcomp).Color = cmapsetE{condcomp}(featI, :);
                bl(condcomp).LineWidth = 1.5;
                bp.FaceColor = cmapsetE{condcomp}(featI, :);
                bp.FaceAlpha = 0.3;
                hold on;
            end

            % Add reference line
            line([-1.500 3.000], repmat(offset, 2), 'Color', cmapsetE{1}(featI, :), 'LineWidth', 1, 'LineStyle', ':');
            axis square;

            % Perform permutation test for statistical significance
            dependent_samples = true;
            p_threshold = 0.1;
            num_permutations = 2000;
            two_sided = true;
            [clusters, p_values, t_sums, permutation_distribution] = permutest(XComp{1, comp}, XComp{2, comp}, dependent_samples, p_threshold, num_permutations, two_sided, 2);

            % Plot significant clusters
            for i = 1:length(clusters)
                if p_values(i) < 0.05
                    line(freqBL.time(clusters{i}), repmat(5 + offset, size(clusters{i})), 'Color', cmapsetE{1}(featI, :), 'LineWidth', 4);
                    text(freqBL.time(clusters{i}(1)), repmat(5 + (offset * 1.2), 1), sprintf('%.3f', p_values(i)), 'Color', cmapsetE{1}(featI, :), 'FontSize', 12);
                end
            end
            savea = [savea, min(p_values, 2)];

            % Set axis limits based on part and component
            if comp ~= 2
                if part == 2
                    axis([-1.500 3.000 -50 5]);
                elseif part == 3
                    axis([-1.500 3.000 -40 20]);
                elseif part == 4
                    axis([-1.500 3.000 -30 10]);
                elseif part == 5
                    axis([-1500 3.000 -30 20]);
                end
            else
                if part == 2
                    axis([-1.500 3.000 -20 20]);
                elseif part == 3
                    axis([-1.500 3.000 -30 20]);
                elseif part == 4
                    axis([-1.500 3.000 -20 20]);
                elseif part == 5
                    axis([-1.500 3000 -30 10]);
                end
            end
            xlim(xlimz{part});
            axis on;
        end
    end
end
set(gcf, 'Position', [546, 97, 1295, 839]);
end