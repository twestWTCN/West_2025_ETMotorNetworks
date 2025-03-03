clear; close all
restoredefaultpath
addpath(genpath(cd))
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                REACHING TO UNDERSTAND CORTICAL NETWORK            %%%%
%%%%                     CHANGES IN ESSENTIAL TREMOR                   %%%%
%%%%                         T. WEST, H. CAGNAN                        %%%%
%%%%                    IMPERIAL COLLEGE LONDON (2025)                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup Path and Configurations
R = projectAddPaths();
R.path.expname = 'EEG_DP_Tremor_V2';
R.path.taskname = {'reaching_task'};
R = setupAnalysisConfig(R);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                             DATA IMPORT                           %%%%  
%%%%                             PREANALYSIS                           %%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Data
% These switches are needed to flick between modalities/sites
R.import.type = 'EEG'; 
R.import.site = 'DL';
R.import.type = 'OPM';
R.import.site = 'UCL';

R.import.subsel = R.import.subselOPAll;
switch R.import.site
    case 'UCL'
        importExperimentalData_OPM(R)
    case 'DL'
        importExperimentalData_EEG(R)
end

%% Preprocess Data
preprocessContinuousData(R);

%% Epoch Data
epochData(R) 

%% Artefact Rejection
R.artRej.subsel = R.import.subselOPAll;
R.artRej.blockpart = {'Rest','Posture','Task'};
R.artRej.meth = 'auto'; % can be 'auto' or 'visual'
artefactRejectData(R)

%% Source Reconstruction
R.locSrcs.subsel = R.import.subselAll;
R.locSrcs.blockpart = {'Rest','Posture','Task'};
% Prepare the leadfields
prepareMRI_EEG_SPM_normalizedGrid(R) % EEG specific
prepareMRI_OPM_normalizedGrid(R) % OPM specific
% Compute the common filters
localizeSourcesCommonFilter_normalizedGrid(R)
% Prepare Vitual Channel Parcellation
computeParcellatedSource_Task(R); % Adjusted to remove OOB channels

%% Spectral Pre-Computation
computeFreqSpectraBank(R,R.import.subselOPAll);

%% Task/trial info
enrichTrialInfo(R,R.import.subselAll)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                         MAIN ANALYSES                            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure(2)  - Task Behaviour
% Perform Task specific Performance
compareControlsVsET_Behaviour(R,R.import.subselAll); 

%% Figure(3) - Tremor Comparisons
compareControlsvsET_Tremor(R,R.import.subselAll);

%% Figure(4) - DICS for Tremor
% Overlap plots
computePostureSimpleTremorDICS_Overlap(R,0)
% Now do contrasts
computePostureSimpleTremorDICS_ContVsET(R,0)
% Now do auxillary
computePostureSimpleTremorDICS_RegressVC(R,0)

%% Figure (5) - DICS for Movement Power
computeTaskSourceContrasts(R,1,1) % EEGSel PatSel

%% Figure (6) - Time Frequency of Different Regions of Interest
compareTaskSpectrograms_VC(R,R.import.subselOPAll); % (Y)

%% Figure (9) - Differences in Beta ET vs Control
compareBetaDesyncPowerVC(R,R.import.subselAll); % (Y)

%% Figure (10) - SPATIOSPECTRAL (PCA)
construct_PCA_matrixSpectrogram(R,R.import.subselAll)
analysePCAComponents(R)
projectPCA2Trials(R,R.import.subselAll)
PCA_RegressComponentsKinematics(R)
