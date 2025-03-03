function compareControlsvsET_Tremor(R, sublist)
% patientVsControlComparison: Analyze tremor properties of patients versus control subjects
% under rest and posture conditions. The function generates summary plots for different
% conditions and performs statistical analyses on the resulting data.

close all;
titnames = {'Rest', 'Posture', 'Cue', 'Reach'};
subI = 0;
patchcols = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]}; % Define colors for plots

% Iterate over each subject in the provided list
for sub = sublist
    subI = subI + 1;
    
    % Iterate over two tasks: Rest and Posture
    for block = 1:2
        retdat.sub = sub;
        retdat.block = block; % Task block
        retdat.cond = []; % Condition (only works for task data)
        retdat.part = []; % Movement part
        retdat.ft = []; % Transformation type
        retdat.ftflag = 1; % FT concatenation flag
        retdat.fileapend = "'pp_arV3'";
        retdat.subfold = 'epoched';
        
        % Retrieve data for analysis
        ftdata_cat = retrieveData(R, retdat);

        % Compute tremor properties
        [peakPow, peakFrq, tremPow, tremVar, tmpPxx] = computeTremorProperties(ftdata_cat, []);
        fxOrgSave(:, block, subI) = squeeze(mean(tmpPxx.powspctrm));
        hz_org = tmpPxx.freq;
        tremPow = log10(tremPow); % Log transform power for analysis

        % Plot histogram of peak frequency
        figure(100);
        if strcmp(sub{1}, 'JM10P031022')
            subplot(1, 2, 1);
            H = histogram(peakFrq, 2:0.5:9, 'Normalization', 'probability');
            H.FaceColor = patchcols{block};
            H.FaceAlpha = 0.6;
            hold on;
            axis square;
            xlabel('Peak Frequency (Hz)');
            ylabel('Probability');
        elseif strcmp(sub{1}, 'TF9P080822')
            subplot(1, 2, 2);
            H = histogram(peakFrq, 2:0.5:9, 'Normalization', 'probability');
            H.FaceColor = patchcols{block};
            H.FaceAlpha = 0.6;
            hold on;
            axis square;
            xlabel('Peak Frequency (Hz)');
            ylabel('Probability');
        end

        % Retrieve subject-specific information
        subTab = retrieveSubjectInfoTab(R, sub);
        TETRAS = table2array(subTab(:, 'TETRAS'));
        if isempty(TETRAS)
            TETRAS = nan;
            limbTREM = nan;
        else
            limbTREM = max(table2array(subTab(:, {'TREM_L', 'TREM_R'})));
        end

        % Calculate summary statistics for each subject
        peakPowX = mean(tremPow);
        peakFrqX = mean(peakFrq);
        powCOV = std(tremPow) / peakPowX;
        frqCOV = std(peakFrq) / peakFrqX;

        if strcmp(sub, 'HH3P160522')
            blockRow(:, block, subI) = nan(1, 6);
            ET(subI) = 2; % Non-ET identifier
        else
            blockRow(:, block, subI) = [peakFrqX, peakPowX, TETRAS, limbTREM, powCOV, frqCOV];
            ET(subI) = strcmp(subTab.Condition{1}, 'ET'); % Assign ET status
        end
    end
end

%% Main Plot - Tremor Analysis for Posture vs Rest
% Main plot to demonstrate tremor differences between ET patients and controls
figure;
close all;
for ip = 1:2
    if ip == 1
        subind = find(ET == 0); % Controls
        cmap = brewermap(24, '*blues');
    else
        subind = find(ET == 1); % ET Patients
        cmap = brewermap(24, '*reds');
    end

    % Plot violin plots for rest and posture conditions
    subplot(2, 2, ip);
    V = violinplot((squeeze(blockRow(2, :, subind)))', {'Rest', 'Posture'}, 'Width', 0.2, 'ViolinColor', cmap(8, :));
    X = [];
    Y = [];
    for i = 1:2
        X(:, i) = V(i).ScatterPlot.XData;
        Y(:, i) = V(i).ScatterPlot.YData;
    end
    P = plot(X(:,:)', Y(:,:)');
    for i = 1:numel(P)
        P(i).Color = cmap(i, :);
    end
    ylabel('Tremor Amplitude (m/s^2)');
    grid on;
    box off;
    axis square;
    a = gca;
    a.XTick = 1:2;
    a.XTickLabel = {'Rest', 'Posture'};
    a.FontSize = 18;
    ylim([-2 1]);
    xlim([0.5 2.5]);

    % Statistical analysis between conditions
    [h, p, ci, stat] = ttest(diff(log(squeeze(blockRow(2, :, subind)))));
    deltaTrem{ip} = 100 .* (squeeze(blockRow(2, 2, subind)) - squeeze(blockRow(2, 1, subind))) ./ squeeze(blockRow(2, 1, subind));
end

%% Plot Frequency Spectrum for ET Group
subplot(2, 2, 3);
XM = mean(squeeze(fxOrgSave(:, 1, find(ET == 1))), 2);
XS = std(squeeze(fxOrgSave(:, 1, find(ET == 1))), [], 2) ./ sqrt(size(find(ET == 1), 2));
[bp(1), pp(1)] = boundedline(hz_org, XM, XS);
bp(1).Color = cmap(8, :);
bp(1).LineStyle = ':';
bp(1).LineWidth = 1.5;
pp(1).FaceColor = cmap(8, :);
pp(1).FaceAlpha = 0.3;

XM = mean(squeeze(fxOrgSave(:, 2, find(ET == 1))), 2);
XS = std(squeeze(fxOrgSave(:, 2, find(ET == 1))), [], 2) ./ sqrt(size(find(ET == 1), 2));
[bp(2), pp(2)] = boundedline(hz_org, XM, XS);
bp(2).Color = cmap(8, :);
bp(2).LineStyle = '-';
bp(2).LineWidth = 1.5;
pp(2).FaceColor = cmap(8, :);
pp(2).FaceAlpha = 0.3;
xlim([2 14]);
grid on;
box off;
axis square;
xlabel('Frequency (Hz)');
ylabel('Power (m/s^2/sqrt(Hz)');
leg = legend(bp, {'Rest', 'Posture'});
leg.Box = 'off';

%% Correlation Analysis between TETRAS Score and Tremor Amplitude
subplot(2, 2, 4);
L = squeeze(blockRow(2, 2, find(ET == 1)));
T = squeeze(blockRow(4, 2, find(ET == 1)));
cmap = brewermap(24, '*reds');
S(1) = scatter(T, L, 75, cmap(8, :), 'filled');
hold on;
[xCalcOut, yCalcOut, b, r2] = plotlinregress(T, L, 0, [3 7]');
plot(xCalcOut, yCalcOut, 'k--', 'LineWidth', 1.5);
grid on;
box off;
axis square;
xlabel('TETRAS Score');
ylabel('Tremor Amplitude (m/s^2)');
set(gcf, 'Position', [281, 464, 1024, 800]);
[r, p] = corr(T(~isnan(T)), L(~isnan(T)), 'type', 'Spearman');

%% Save Tremor Assessment Data to File
T = table(sublist', double(ET'), squeeze(blockRow(2, 2, :)), squeeze(blockRow(2, 1, :)), vertcat(deltaTrem{:}), squeeze(blockRow(4, 1, :)), squeeze(blockRow(3, 1, :)), squeeze(blockRow(6, 1, :)));
T.Properties.VariableNames = {'Name', 'ET', 'Postural Tremor', 'Rest Tremor', 'Rest2Post', 'TETRAS_Trem', 'TETRAS_Perf', 'PostureTremPeak'};
writetable(T, 'TremorAssessmentData.xlsx');
T = retrieveTremorInfoTab(R, sub);
end