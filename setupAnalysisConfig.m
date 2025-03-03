function R = setupAnalysisConfig(R)
%% This function sets up the analysis configurations and parameters

%% PRESERVE THESE LISTS
R.import.subselAll = {
    'AK10C170822','DP8C300522','GF4C050422',...
    'RK9C150622','SF5C060422',...
    'SS12C120122','TF11C071022','TW7C250522',...
    'HA13C130223','VR14C150323','MM15C170323' ...% Controls first
    ... % Now ETs
    'BP8P120722','BS13P101022','HB5P100622',...
    'HH3P160522','JM10P031022','JS6P140622',...
    'MBP7P170622','SK11P051022','TF9P080822',...
    'HW14P141022','JK16P270223','DM15P301122'...
    };

% Controls
R.import.subselCont = {'AK10C170822','DP8C300522','GF4C050422',...
    'RK9C150622','SF5C060422',...
    'SS12C120122','TF11C071022','TW7C250522',...
    'HA13C130223','VR14C150323','MM15C170323'};
% Patients
R.import.subselPat = {'BP8P120722','BS13P101022','HB5P100622',...
    'HH3P160522','JM10P031022','JS6P140622',...
    'MBP7P170622','SK11P051022','TF9P080822',...
    'HW14P141022','JK16P270223','DM15P301122'};

% EEG Controls
R.import.subselOPAll = {'OP00043','OP00055','OP00067','OP00069','OP00068','OP00103','OP00134','OP00137','OP00139'};
R.import.subselOPPat  = {'OP00068','OP00103','OP00134','OP00139'};
R.import.subselOPCont = {'OP00043','OP00055','OP00067','OP00069','OP00137'};

% OPM Specifics
R.opmsense.trigchan = {'NI-TRIG-1'};
R.opmsense.ACCchan = {'NI-TRIG-2','NI-TRIG-3','NI-TRIG-4'};

%% Import data
% This will import raw data files, convert to fieldtrip, and then split up
% data into the respective blocks: rest,posture, and task using minimal
% user input.
R.import.condsel = 1:4;
R.import.repsel = [{1:4}]; %These are the recording repetitions; This MUST match subsel codes {1:2}
R.import.blockpart = {'Rest','Posture','Task'};
R.import.condname = {'High UC, Low Prec','Low UC, Low Prec','High UC, High Prec','Low UC, High Prec'};
R.import.condname_short = {'HUC,LRG','LUC,LRG','HUC,SML','LUC,SML'};
R.import.condname_groups{1} = {'HUC','LUC'};
R.import.condname_groups{2} = {'LRG','SML'};
R.import.rs_fs = 512; % resample speed;
R.import.repsel = {[]}; % Full
R.import.TAK.visplit = [0 0 0 0];

R.opmsense.trigchan = {'NI-TRIG-1'};
R.opmsense.ACCchan = {'NI-TRIG-2','NI-TRIG-3','NI-TRIG-5'};

%% Preprocess
R.import.blockpart ={'Rest','Posture','Task'};
R.preprocess.subsel= R.import.subselAll;
R.preprocess.ICA.flag = 1; % compute ICA components fresh

%% Epoch
R.epoch.blockpart = {'Rest','Posture','Task'};
R.epoch.subsel = R.import.subselAll;
R.epoch.accBadFlag = {zeros(size(R.epoch.subsel))};
R.epoch.accBadFlag{1}(strcmp(R.epoch.subsel,'HH3P160522')) = 1;
R.epoch.length = 2; % Epoch length for steady state (s)
R.epoch.burn = 1; %Number of Epochs to remove either end (N trials)
R.epoch.names = {'TaskRest','Posture','Prep','Exec','Hold','FullSeq'};
R.epoch.blockpart = {'Rest','Posture','Task'};
R.epoch.manfresh = 0;
% Trial length definition
% trial 1: search for movement -500 to 2500ms posture;
% trial 2: define as 0 to 1500ms arrow cue
% trial 3: search for movement 0 to 2500ms go cue
% trial 4: define as 0 to 2500ms collision mark
% R.epoch.trans = {[-1500 4000],[-1500 4000],[-1500 4000],[-1500 4000],[-3000 4000],[0 0]}; % defines actual length of window
R.epoch.trans = {[-3000 4000],[-3000 4000],[-3000 4000],[-3000 4000],[-3000 4000],[0 0]}; % V3 - allows better spectrogram average
R.epoch.block = {[0 1500],[0 1500],[0 1500],[0 1500],[0 1500],[0 0]}; % defines actual length of window
R.epoch.searchwin = {[0 2500],...
    [0 2500],... % for postural hold, allow a prepotent response
    [0 0],...,
    [-1500 6000],...
    [-2000 3000],...
    [0 0]}; % Search window (ms for locking to movement)

%% Artefact Rejection
R.artRej.blockpart = {'Rest','Posture','Task'};
R.artRej.meth = 'auto'; % can be 'auto' or 'visual'
% % R.artRej.maxZScoreEps = 15; % Max Z-Score threshold
% % R.artRej.maxZScoreEpsFlat = 1.5;
% % R.artRej.maxZScoreEpsCh = 12;
% % R.artRej.maxZScoreMusc = 18;

R.artRej.maxZScoreEps = 10; % Max Z-Score threshold
R.artRej.maxZScoreEpsFlat = 1.5;
R.artRej.maxZScoreEpsCh = 10;
R.artRej.maxZScoreMusc = 10;

%% Correct Movement Axes
R.corrMoveAxes.blockpart = {'Task'};
R.corrMoveAxes.subsel = R.import.subselAll;

%% PLOTTING
% Setup matlab defaults plotting
set(0,'defaultAxesFontSize',18)
set(0,'DefaultAxesFontName','Arial')
R.plot.cmap.lplots = linspecer(6);

%% Export options
R.export.pathonly = 0;