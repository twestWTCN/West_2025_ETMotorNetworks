function compareControlsVsET_Behaviour(R,sublist)
% analyseMovementET_vsControl: Analyze movement kinematics and accelerometer data
% for Essential Tremor (ET) versus control subjects under various conditions.
% This function processes data from multiple subjects and saves analysis results.

close all
fresh = 1;

if fresh
    % Initialize main variables
    TAB = []; TABZs = []; subI = 0; magSave = []; accSave = [];  magSpaceSave = []; accSpaceSave = [];

    % Iterate over each subject in the provided list
    for sub = sublist
        subI = subI + 1;
        TAB_sub = [];

        % Iterate over four conditions for each subject
        for cond = 1:4

            %% Analyse the Hold phase
            retdat.sub = sub;
            retdat.block = 3; % task type
            retdat.cond = cond; % condition number
            retdat.part = 5; % movement part
            retdat.ft = 2; % transformation type
            retdat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']"; % Use non-artefact rejected version
            retdat.subfold = "['Condition' num2str(cond) '_epoched']";
            retdat.ftflag = 1; % nonFT concatenation flag
            % Retrieve data for analysis
            ftdata_cat = retrieveData(R, retdat);

            % Identify accelerometer and position signal channels
            accSig = strcmp(ftdata_cat.hdr.chantype, 'Accelerometer');
            motiveSig = strcmp(ftdata_cat.hdr.chantype, 'MotivePos');
            fsamp = ftdata_cat.fsample;

            % Double-check and adjust position signal labels if necessary
            labs = ftdata_cat.label(motiveSig);
            if labs{1}(1) ~= 'x'
                motiveSig = strcmp(ftdata_cat.label, 'xFinger1Pos') | strcmp(ftdata_cat.label, 'yFinger1Pos') | strcmp(ftdata_cat.label, 'zFinger1Pos');
                if isempty(motiveSig)
                    motiveSig = strcmp(ftdata_cat.label, 'PosX') | strcmp(ftdata_cat.label, 'PosY') | strcmp(ftdata_cat.label, 'PosZ');
                end
            end

            % Initialize variables for velocity and acceleration calculations
            magnitude_velocityNorm = [];
            magnitude_accNorm = [];
            magnitude_velocitySpace = nan(size(ftdata_cat.trial, 2), 4000);
            magnitude_accSpace = nan(size(ftdata_cat.trial, 2), 4000);

            % Loop through trials and compute kinematics
            for trn = 1:size(ftdata_cat.trial, 2)
                motiveData = ftdata_cat.trial{trn}(motiveSig, :);
                accData = signal2accel(ftdata_cat.trial{trn}(accSig, :));
                magnitude_acc = sqrt(sum(accData.^2, 1));

                timeInd = find(ftdata_cat.time{trn} >= -3, 1):find(ftdata_cat.time{trn} > 2, 1);
                XYZ = motiveData(:, timeInd);
                [~, ~, ~, ~, ~, ~, magnitude_velocity, moveindex] = computeKinematics(XYZ, fsamp, 0.2, 0.3);

                % Resample to standard size
                Xbase = 1:numel(magnitude_velocity);
                baseVector = linspace(Xbase(1), Xbase(end), 2048);
                magnitude_velocitySpace(trn, 1:numel(magnitude_velocity)) = magnitude_velocity;
                X = magnitude_acc(timeInd(moveindex)); X = X - X(1);
                magnitude_accSpace(trn, 1:numel(X)) = X;
                magnitude_velocityNorm(:, trn) = interp1(Xbase, magnitude_velocity, baseVector);
                magnitude_accNorm(:, trn) = interp1(Xbase, X, baseVector);
            end

            % Save calculated data
            magSave = [magSave magnitude_velocityNorm];
            accSave = [accSave magnitude_accNorm];
            magSpaceSave = [magSpaceSave; magnitude_velocitySpace];
            accSpaceSave = [accSpaceSave; magnitude_accSpace];

            magVelSave(:, cond, subI) = mean(magnitude_velocityNorm, 2);
            magaccSave(:, cond, subI) = mean(magnitude_accNorm, 2);
            TAB_sub = [TAB_sub; ftdata_cat.trialinfo];
        end

        % Append results for each subject
        TAB = [TAB; TAB_sub]; % Raw data
        TABZs = [TABZs; nanzscore(TAB_sub, 1)]; % Z-scored data per subject
    end

    % Define variable names for output tables
    parList = {'Cond', 'Dir', 'Success', 'ExecTime', 'HoldTime', 'ReachQ', 'SubI', 'ET', ...
        'CorrectTar', 'errorTar', 'holdError', 'holdVar', 'goCueOnset', 'holdDur', ...
        'tremPeakPow', 'tremFreq', 'tremPow', 'tremVar', 'mean_velocity', 'max_velocity', ...
        'pathlength', 'mov_dur', 'smooth_index', 'trialDeviance', 'holdStability', 'reachDur', ...
        'reachRT', 'holdRT'};

    % Create tables and assign column names
    matTAB = array2table(double(TAB));
    matTAB.Properties.VariableNames = parList;
    matTAB.ReachDur = matTAB.HoldTime - matTAB.ExecTime;
    matTABzs = array2table(double(TABZs));
    matTABzs.Properties.VariableNames = parList;

    % Convert relevant columns to categorical data types
    matTAB.Cond = categorical(matTAB.Cond);
    matTAB.SubI = categorical(matTAB.SubI);
    matTAB.UC = categorical(matTAB.Cond == '1' | matTAB.Cond == '3');
    matTAB.Prec = categorical(matTAB.Cond == '1' | matTAB.Cond == '2');
    matTAB.ET = categorical(matTAB.ET);
    matTAB.ReachQ = categorical(matTAB.ReachQ);

    matTABzs.Cond = matTAB.Cond;
    matTABzs.SubI = matTAB.SubI;
    matTABzs.UC = matTAB.UC;
    matTABzs.Prec = matTAB.Prec;
    matTABzs.ET = matTAB.ET;
    matTABzs.ReachQ = matTAB.ReachQ;

    % Save analysis results
    save('analyseMovementGoalTabHold', 'matTAB', 'matTABzs', 'magSave', 'accSave', 'magVelSave', 'magaccSave', 'magSpaceSave', 'accSpaceSave')
else
    % Load previously saved analysis results
    load('analyseMovementGoalTabHold', 'matTAB', 'matTABzs', 'magSave', 'accSave', 'magVelSave', 'magaccSave', 'magSpaceSave', 'accSpaceSave')
end

%% Further Analysis and Visualization
% Parameters for analysis
parListNeat = {'Reach Q', 'Reach RT', 'Mean Vel', 'Max Vel', 'Pathlength', ...
    'Move Duration', 'Trial Deviance', 'Hold Duration', 'Hold Stability', 'Hold RT','Tremor Power'};
parList = {'ReachQ', 'reachRT', 'mean_velocity', 'max_velocity', 'pathlength', ...
    'mov_dur', 'trialDeviance', 'holdDur', 'holdStability', 'holdRT', 'tremPeakPow'};

% Compute subject-specific statistics
for parI = 1:numel(parList)
    for cond = 1:4
        for sub = unique(matTAB.SubI)'
            % Determine condition indices
            if cond == 1
                condind = matTAB.Prec == 'true';
            elseif cond == 2
                condind = matTAB.Prec == 'false';
            elseif cond == 3
                condind = matTAB.UC == 'true';
            elseif cond == 4
                condind = matTAB.UC == 'false';
            end

            % Determine subject index
            subind = matTAB.SubI == sub;
            combind = subind & condind;

            % Calculate mean parameter statistics
            if parI < 11
                partStat{parI}(cond, sub) = nanmean(double(matTAB.(parList{parI})(combind)));
            else
                partStat{parI}(cond, sub) = nanmean(log10(matTAB.(parList{parI})(combind))); % log scale the tremor parameter
            end

            % Save magnitude velocity and acceleration data
            magVelSpaceG(:, cond, sub) = nanmean(magSpaceSave(combind, :), 1);
            magVelG(:, cond, sub) = mean(magSave(:, combind), 2);
            ETG(cond, sub) = double(matTAB.ET(find(combind, 1)));
            CondG(cond, sub) = cond;
            SubG(cond, sub) = double(sub);
        end
    end
    % Subtract subject average
    partStat{parI} = partStat{parI} - mean(partStat{parI}, 1);
end

% Partition by parameters for further analysis
for parI = 1:numel(parList)
    for sub = unique(matTAB.SubI)'
        for partit = 1:3
            % Find subject indices
            subind = matTAB.SubI == sub;
            subind = find(subind);
            X = double(matTAB.(parList{parI})(subind));

            % Partition trials based on parameter values
            if partit == 1
                ls = find(X <= ceil(nanmedian(X)));
            elseif partit == 2
                ls = find(X > nanmedian(X));
            elseif partit == 3
                ls = 1:numel(X);
            end

            % Calculate mean velocity and acceleration profiles for each partition
            magVelSpaceG_part(:, partit, parI, sub) = nanmean(magSpaceSave(subind(ls), :), 1);
            magVelG_part(:, partit, parI, sub) = mean(magSave(:, subind(ls)), 2);
            accVelSpaceG_part(:, partit, parI, sub) = nanmean(abs(accSpaceSave(subind(ls), :)), 1);
            accVelG_part(:, partit, parI, sub) = mean(abs(accSave(:, subind(ls))), 2);
        end
    end
end

%% Plot the Velocity Profiles: Simple ET vs Control
cmap = brewermap(2, 'Set1');
ytit = {'Velocity (M/s)', 'Acceleration (M/s^2)'};

for i = 1:2
    subplot(1, 2, i)
    for ET = 1:2
        % Select subjects based on ET group
        if ET == 1
            subsel = 1:11;
        else
            subsel = 12:23;
        end

        % Select appropriate data for velocity or acceleration
        if i == 1
            Z = magVelG_part(:, :, 1, subsel);
        else
            Z = accVelG_part(:, :, 1, subsel);
        end

        % Calculate mean and standard error across subjects
        X = squeeze(nanmean(nanmean(Z, 3), 4));
        XS = squeeze(nanstd(Z, [], 4)) ./ sqrt(numel(subsel));
        ccmap = cmap(ET, :);

        % Plot bounded line for partition 3
        for partit = 3
            fsamp = 512;
            tvec = linspace(0, 100, size(X, 1));
            [l(partit, ET), b(partit, ET)] = boundedline(tvec, X(:, partit), XS(:, partit));
            l(partit, ET).Color = ccmap;
            b(partit, ET).FaceColor = ccmap;
            b(partit, ET).FaceAlpha = 0.2;
            hold on
        end
        hold on
        axis square;
        box off;
        grid on;
    end
    ylabel(ytit{i})
    xlabel('Normalized time (%)')
    legend(l(3, :), {'Control', 'ET'}, 'box', 'off')
end
%% Box Plots Analysis
% Plotting parameter statistics across conditions
R.import.condname_short = {'LRG', 'SML', 'HUC', 'LUC'};
axlims = {[-0.15 0.15], [-500 500], [-0.05 0.05], [-0.1 0.1], [-2.5 2.5], [-150 150], [-0.05 0.05], [-600 600], 1e-4.*[-1 1], [-500 500]};
figure(1231); clf;
parSel = [2, 3, 5, 10];

for parI = 1:numel(parSel)
    subplot(1, 4, parI);
    [B, S, p] = scatterBoxChart(R.import.condname_short, partStat{parSel(parI)}(:), CondG(:), double(ETG(:)), cmap);

    a = gca;
    a.XTickLabel = R.import.condname_short;
    ylabel(parListNeat{parSel(parI)});
    axis square; box off; a.FontSize = 12;
    set(gcf, 'Position', [133, 392, 1400, 600]);

    % Create a table for each parameter and store it
    TAB = array2table([partStat{parSel(parI)}(:), SubG(:), CondG(:) == 1, CondG(:) == 3, ETG(:), partStat{11}(:)], ...
        'VariableNames', {parList{parSel(parI)}, 'SubI', 'Prec', 'UC', 'ET', 'tremPeakPow2'});
    partable{parI} = TAB;
end
