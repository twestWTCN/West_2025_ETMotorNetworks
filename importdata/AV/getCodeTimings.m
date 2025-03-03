function [decoder, SMRtrialdef] = getCodeTimings(CVM, ampCoder, trialdetails, splitType, blockSize, minbreak)
    % CODETIMINGSEEG: Decodes the stages of a trial from coded signals, splits the data
    % into relevant blocks, and returns timing definitions for further analysis.
    % Parameters:
    %   CVM - Calibration values from trigger calibration measurement
    %   ampCoder - Coder signal from amplifier for movement/calibration trial
    %   trialdetails - Details about the trial to be appended
    %   splitType - Method of trial splitting ('blockSplit', 'taskSplit', 'condSplit')
    %   blockSize - Minimum block sizes for different stages
    %   minbreak - Minimum allowed break length between trials

    %% Step 1: Decode Trial Stages from Amplifier Data
    ampCoder(isnan(ampCoder)) = 0; % Replace NaNs with zeros
    ampDeCoder = abs(ampCoder - CVM); % Calculate absolute difference
    [dError, decoder] = min(ampDeCoder, [], 1); % Decode to closest calibration value
    decoder(dError > 0.01) = nan; % Filter out values above error threshold

    %% Step 2: Generate Trial Definitions (trialdef)
    for trialstage = 1:numel(CVM)
        % Identify timings at which each trial stage exists
        X = find(decoder == trialstage);
        blocks = SplitVec(X, 'consecutive');
        trialdef{trialstage} = zeros(2, numel(blocks));

        % Loop through repetitions to identify block boundaries
        nmx = []; % Track block sizes
        for SR = 1:numel(blocks)
            nmx(SR) = numel(blocks{SR});
            trialdef{trialstage}(:, SR) = [blocks{SR}(1); blocks{SR}(end)]; % Start and end samples of each block
        end

        %% Step 2.1: Remove Small Breaks Between Blocks
        trBreak = [inf, diff(trialdef{trialstage}(2, :))];
        p = 0; % Counter to adjust for changes in length
        for i = find(trBreak < minbreak(trialstage))
            p = p + 1;
            A = trialdef{trialstage}(1, :);
            B = trialdef{trialstage}(2, :);
            A(i - p) = [];
            B(i - p - 1) = [];
            nmx(i - p) = [];
            trialdef{trialstage} = [A; B];
        end

        %% Step 2.2: Remove Blocks that Are Too Small
        trialdef{trialstage}(:, nmx < blockSize(trialstage)) = [];
        nmx(nmx < blockSize(trialstage)) = [];
    end

    %% Step 3: Fix Initial Null Blocks (If Any)
    if trialdef{1}(1, 1) < trialdef{3}(1, 1)
        trialdef{1}(:, 1) = []; % Remove leading null block
    end

    %% Step 4: Split Trials Based on splitType
    switch splitType
        case 'blockSplit'
            % Block Split into Rest, Posture, and Task
            naming = {'Rest', 'Posture', 'Task'};
            restdef = trialdef{2}(:, 1);
            postdef = trialdef{3}(:, 1);
            taskdef = [trialdef{2}(2, 1), trialdef{1}(1, 1);
                       trialdef{1}(2, 1), trialdef{1}(1, 2);
                       trialdef{1}(2, 2), trialdef{1}(1, 3);
                       trialdef{1}(2, 3), trialdef{4}(end)];

            trialdef = {restdef, postdef, taskdef};
            SMRtrialdef = [naming; trialdef]';
            SMRtrialdef = [SMRtrialdef, repmat(trialdetails, size(trialdef, 2), 1)];

        case 'taskSplit'
            % Split into Motor Rest, Preparation, Execution, and Hold
            naming = {'MotorRest', 'MotorPrep', 'MotorExec', 'Hold'};
            motorRest = trialdef{1};
            motorPrep = trialdef{2};
            motorExec = trialdef{3};
            motorHold = trialdef{4};

            trialdef = {motorRest, motorPrep, motorExec, motorHold};
            SMRtrialdef = [naming; trialdef]';
            SMRtrialdef = [SMRtrialdef, repmat(trialdetails, size(trialdef, 2), 1)];

        case 'condSplit'
            % Split by Conditions ('1', '2', '3', '4')
            naming = {'1', '2', '3', '4'};
            catCond = [trialdef{:}];
            condDef = cell(1, numel(naming));

            % Loop through conditions and assign trials based on proximity
            for cond = 1:numel(naming)
                TD = trialdef{cond}(1, :);
                tmpDef = [];
                for p = 1:numel(TD)
                    trDist = TD(p) - catCond(1, :); % Find distance to next trial
                    trDist(trDist >= 0) = nan;
                    if any(trDist)
                        [~, k] = max(trDist);
                        tmpDef = [tmpDef, [TD(p); catCond(1, k)]];
                    else
                        tmpDef = [tmpDef, [TD(p); size(ampCoder, 2)]];
                    end
                end
                condDef{cond} = tmpDef;
            end

            trialdef = condDef;
            SMRtrialdef = [naming; trialdef]';
            SMRtrialdef = [SMRtrialdef, repmat(trialdetails, size(trialdef, 2), 1)];
    end
end
