function AVData = importEngAVData(filename,fsamp)
% Script for importing AV Data (from Unity program) for UK Computers

opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Condition", "Coder", "Timer", "xScreenPosition", "yScreenPosition", "trackingX", "trackingY", "trackingZ"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
AVData = readtable(filename, opts);

time_AV = AVData.Timer;
data_AV = [AVData.Condition AVData.Coder AVData.xScreenPosition AVData.yScreenPosition AVData.trackingX AVData.trackingY AVData.trackingZ];

% Remove incomplete rows
data_AV(isnan(time_AV),:) = [];
time_AV(isnan(time_AV),:) = [];

time_AV(any(isnan(data_AV),2),:) = [];
data_AV(any(isnan(data_AV),2),:) = [];

% Resample to Target
tartime_AV = time_AV(1):1./fsamp:time_AV(end);
% Ensure no duplicate points
[C,IA] = unique(time_AV);
time_AV = time_AV(IA);
data_AV = data_AV(IA,:);

trials = interp1(time_AV,data_AV,tartime_AV');
time_AV = tartime_AV;

%% Setup FT data structure
AVData = [];
AVData.fsample = fsamp;
AVData.trial{1} = trials';
AVData.time{1} = time_AV;
AVData.label = {'Condition','Coder','xScreen','yScreen','xFinger1Pos','yFinger1Pos','zFinger1Pos'};
AVData.chantype = {'Info','Trigger','UnityPos','UnityPos','MotivePos','MotivePos','MotivePos'};

save([fileparts(filename) '\AVData.mat'],'AVData')
end