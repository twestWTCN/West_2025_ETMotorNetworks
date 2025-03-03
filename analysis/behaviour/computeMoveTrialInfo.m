function ftdata_cat = computeMoveTrialInfo(ftdata_cat, vargin)
% Compute movement trial information and update FieldTrip data structure
% Inputs:
%   ftdata_cat - FieldTrip data structure
%   vargin - Variable input arguments
%      R - structure containing paths and experiment settings
%      sub - subject identifier
%      subI - subject index
%      plotOp - plotting option flag
%      part - task part identifier

% Extract variable inputs
R = vargin{1};
sub = vargin{2};
subI = vargin{3};
plotOp = vargin{4};
part = vargin{5};

% Load the model trajectories
fileappend = [R.epoch.names{part} '_modelTraj'];
if part == 2 || part == 4 || part == 5
    modelTraj = loadExpData(R, sub{1}, 'Task', [], fileappend, ['AllConditions']);
end

% Key codes for events:
%   1: pretrial rest               = 0.2v
%   2: pretrial posture            = 0.5v
%   3: movement cue (black arrows) = 1.5v
%   4: go signal (green arrows)    = 2.0v
%   5: hold (collision)            = 2.5v
%   6: successful hold (pop)       = 3.0v

% Actual centre of mass for conditions
pos{1} = [...
    1.4142         0   -1.4142   -2.0000   -1.4142         0    1.4142    2.0000;
    1.5979    2.1837    1.5979    0.1837   -1.2304   -1.8162   -1.2304    0.1837];

pos{2} = [...
    1.4142         0   -1.4142   -2.0000   -1.4142         0    1.4142    2.0000;
    1.5061    2.0919    1.5061    0.0919   -1.3223   -1.9081   -1.3223    0.0919];

% Identify condition from trial info
cond = unique(ftdata_cat.trialinfo(:, 1));
if numel(cond) > 2
    error('This data has mixed conditions, compute movement trial info won''t work!');
end

% Set target center based on condition
if cond <= 2
    unityTarPos = pos{1};
else
    unityTarPos = pos{2};
end

% Get indices of movement signals
srcs = {'MotivePos', 'UnityPos', 'unitypos', 'Accelerometer', 'Coder'};
unitySig = strcmp(ftdata_cat.hdr.chantype, srcs{2}) | strcmp(ftdata_cat.hdr.chantype, srcs{3});
motiveSig = strcmp(ftdata_cat.hdr.chantype, srcs{1});
accSig = strcmp(ftdata_cat.hdr.chantype, srcs{4});
coderSig = strcmp(ftdata_cat.label, 'Coder');
fsamp = ftdata_cat.fsample;

% Double-check labels
labs = ftdata_cat.label(motiveSig);
if labs{1}(1) ~= 'x'
    motiveSig = strcmp(ftdata_cat.label, 'xFinger1Pos') | strcmp(ftdata_cat.label, 'yFinger1Pos') | strcmp(ftdata_cat.label, 'zFinger1Pos');
    if isempty(motiveSig)
        motiveSig = strcmp(ftdata_cat.label, 'PosX') | strcmp(ftdata_cat.label, 'PosY') | strcmp(ftdata_cat.label, 'PosZ');
    end
    unitySig = strcmp(ftdata_cat.label, 'xScreen') | strcmp(ftdata_cat.label, 'yScreen');
end

% Initialize trial-related variables
actualTarget = []; errorTarget = []; holdError = []; holdVar = [];
holdDur = []; goCueOnset = []; tremPostPow = [];
peakPow = []; peakFrq = []; tremPow = []; tremVar = [];
mean_velocity = []; max_velocity = []; pathlength = [];
mov_dur = []; smooth_index = []; holdStability = [];
reachDur = []; reachRT = []; holdRT = []; trialDeviance = [];

% Loop through trials
for trn = 1:numel(ftdata_cat.trial)
    % Common variables for the current trial
    timeInd = find(ftdata_cat.time{trn} >= 0, 1):find(ftdata_cat.time{trn} > 2, 1);
    timeData = ftdata_cat.time{trn};
    coderData = ftdata_cat.trial{trn}(coderSig, :);
    motiveData = ftdata_cat.trial{trn}(motiveSig, :);
    unityData = ftdata_cat.trial{trn}(unitySig, :);
    accData = ftdata_cat.trial{trn}(accSig, :);

    % Set default values to NaN
    goCueOnset(trn) = nan;
    holdDur(trn) = nan;
    actualTarget(trn) = nan;
    errorTarget(trn) = nan;
    holdError(trn) = nan;
    holdVar(trn) = nan;
    trialDeviance(trn) = nan;
    mean_velocity(trn) = nan;
    max_velocity(trn) = nan;
    pathlength(trn) = nan;
    mov_dur(trn) = nan;
    smooth_index(trn) = nan;
    holdStability(trn) = nan;
    reachDur(trn) = nan;
    reachRT(trn) = nan;
    holdRT(trn) = nan;
    peakPow(trn) = nan;
    peakFrq(trn) = nan;
    tremPow(trn) = nan;
    tremVar(trn) = nan;
    tremPostPow(trn) = nan;

    baseVectorN = 2048;
    if part == 2
        % Posture Kinematics
        timeInd = find(ftdata_cat.time{trn} >= -0.5, 1):find(ftdata_cat.time{trn} > 2, 1);
        XYZ = motiveData(:, timeInd);
        [mean_velocity(trn), max_velocity(trn), pathlength(trn), mov_dur(trn), smooth_index(trn), holdStability(trn), magnitude_velocity, moveIndex] = computeKinematics(XYZ, fsamp, 0.2, 0.3);

        % Project the magnitude velocity to the common base vector
        Xbase = 1:numel(magnitude_velocity);
        baseVector = linspace(Xbase(1), Xbase(end), baseVectorN);
        accTraj = interp1(Xbase, magnitude_velocity, baseVector);
        trialDeviance(trn) = rmse(accTraj, modelTraj', 2);
        reachDur(trn) = 1000 * length(moveIndex) / fsamp;

        % Compute tremor properties
        latency = [ftdata_cat.time{trn}(timeInd(moveIndex(end))) ftdata_cat.time{trn}(end) trn];
        if diff(latency(1:2)) > 1
            [peakPow(trn), peakFrq(trn), tremPow(trn), tremVar(trn)] = computeTremorProperties(ftdata_cat, latency, 0);
            tremPow(trn) = log10(tremPow(trn));
        end
    elseif part == 3
        % Identify Go Cue
        goCue = findGoCue(R, coderData);
        goCueOnset(trn) = isempty(goCue) * inf + ~isempty(goCue) * ftdata_cat.time{trn}(goCue);

        % Compute kinematics
        timeInd = find(ftdata_cat.time{trn} < 0, 1):find(ftdata_cat.time{trn} < 2, 1, 'last');
        XYZ = motiveData(:, timeInd);
        [mean_velocity(trn), max_velocity(trn), pathlength(trn), mov_dur(trn), smooth_index(trn), holdStability(trn), magnitude_velocity, moveIndex, holdIndex] = computeKinematics(XYZ, fsamp, 0.1, 0.8);

        % Compute tremor properties
        latency = [0 goCueOnset(trn) trn];
        if diff(latency(1:2)) > 1
            [peakPow(trn), peakFrq(trn), tremPow(trn), tremVar(trn)] = computeTremorProperties(ftdata_cat, latency, 0);
            tremPow(trn) = log10(tremPow(trn));
        end
    elseif part == 4 || part == 5
        % Reach Analysis
        goCue = findGoCue(R, coderData);
        goCueOnset(trn) = isempty(goCue) * inf + ~isempty(goCue) * ftdata_cat.time{trn}(goCue);

        % Reach Kinematics
        timeInd = find(ftdata_cat.time{trn} >= -0.5, 1):find(ftdata_cat.time{trn} < 4, 1, 'last');
        XYZ = motiveData(:, timeInd);
        [mean_velocity(trn), max_velocity(trn), pathlength(trn), mov_dur(trn), smooth_index(trn), holdStability(trn), magnitude_velocity, moveIndex, holdIndex] = computeKinematics(XYZ, fsamp, 0.2, 0.3);

        % Calculate additional metrics
        reachRT(trn) = 1000 * (ftdata_cat.time{trn}(timeInd(moveIndex(1))) - goCueOnset(trn));
        reachDur(trn) = 1000 * length(moveIndex) / fsamp;
        Xbase = 1:numel(magnitude_velocity);
        baseVector = linspace(Xbase(1), Xbase(end), baseVectorN);
        accTraj = interp1(Xbase, magnitude_velocity, baseVector);

        % Compute tremor properties for posture
        latency = [-3 -1.5 trn];
        [~, ~, tremPostPow(trn), ~] = computeTremorProperties(ftdata_cat, latency, 0);
        tremPostPow(trn) = log10(tremPostPow(trn));

        % Handle reach success and errors
        if ~isnan(holdIndex)
            targetArrow = ftdata_cat.trialinfo(trn, 2) / 45;
            if isempty(unityData)
                actualTarget(trn) = targetArrow;
                distXY = nan;
                tarInd = nan;
                subflag(subI) = 1;
            else
                [actualTarget(trn), distXY] = getActualTarget(unityData(:, holdIndex), cond);
                [~, tarInd] = nanmin(nanmin(distXY));
            end
            trialDeviance(trn) = rmse(accTraj, modelTraj(:, actualTarget(trn))', 2);
            correctTarget(trn) = targetArrow == tarInd;
            holdDur(trn) = 1000 * length(holdIndex) / fsamp;
            distXY(distXY > 2) = nan;
            holdError(trn) = nanmean(distXY);
            holdVar(trn) = nanstd(distXY);
            holdRT(trn) = 1000 * (ftdata_cat.time{trn}(timeInd(holdIndex(1))) - goCueOnset(trn));

            % Calculate tremor properties during hold
            latency = [ftdata_cat.time{trn}(timeInd(holdIndex(1))), ftdata_cat.time{trn}(timeInd(holdIndex(end))), trn];
            if diff(latency(1:2)) > 1
                [peakPow(trn), peakFrq(trn), tremPow(trn), tremVar(trn)] = computeTremorProperties(ftdata_cat, latency, 0);
                tremPow(trn) = log10(tremPow(trn));
            end
        end
    end
end

% Update trial information in the FieldTrip structure
if sub{1}(end - 6) == 'C'
    ET = 0;
elseif sub{1}(end - 6) == 'P'
    ET = 1;
else
    ET = 3;
end

ftdata_cat.trialinfo = [ftdata_cat.trialinfo(:, 1:6), repmat(subI, numel(holdVar), 1), repmat(ET, numel(holdVar), 1),...
    actualTarget', errorTarget', holdError', holdVar', goCueOnset', holdDur',...
    peakPow', peakFrq', tremPow', tremVar', mean_velocity', max_velocity',...
    pathlength', mov_dur', smooth_index', trialDeviance', holdStability', reachDur',...
    reachRT', holdRT', tremPostPow'];

ftdata_cat.trialinfo = double(ftdata_cat.trialinfo);
ftdata_cat.cfg.trialInfoHeader = {'Cond', 'Dir', 'Success', 'ExecTime', 'HoldTime', 'ReachQ', 'SubI', 'ET',...
    'CorrectTar', 'errorTar', 'holdError', 'holdVar', 'goCueOnset', 'holdDur',...
    'tremPeakPow', 'tremFreq', 'tremPow', 'tremVar', 'mean_velocity', 'max_velocity',...
    'pathlength', 'mov_dur', 'smooth_index', 'trialDeviance', 'holdStability', 'reachDur',...
    'reachRT', 'holdRT', 'postremPow'};

end

function goCue = findGoCue(R, coder)
% Find the go cue based on modality type
switch R.import.type
    case 'EEG'
        goCue = find(((coder - 0.020).^2) < 1e-6, 1, 'first');
    case 'OPM'
        goCue = find(((coder - 2.000).^2) < 1e-6, 1, 'first');
end
end
