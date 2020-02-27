%% Setup Trial Table and Kinematics Index File
load('C:\Users\dsarma\Box Sync\IM EMG for R01\Kinematics_Organized.mat'); %Carl's Kinematics based Movement Indexes

trials = readtable('R:\data_raw\human\emg_decoding\EMGC14\20170424_s15\Exports\metaTable\EMG\EMG_30kHz\emgraw_metaTrial.csv','ReadRowNames',true); %Trials table form HAPTIX
trials(9,:) = []; %Remove Dead File/Trial
Fs = 30e3;

for actionCase = 1:8
    switch actionCase
    case 1 %% Wrist - Pro/Sup - Neutral Position - Flat
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0053.mat'
        muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Wrist.PS.Neut.Flat.Idx_25'-4500;  
        
    case 2 %% Wrist - Flx/Ext - Wrist Pronated - Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0044.mat'
        muscleName = {'EPL', 'APL', 'FPL', 'SUP', 'PTER', 'ED3', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Wrist.FE.Pro.Flat.Idx_25'-4500;
        
    case 3 %% D1 - Flx/Ext - Wrist Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0013.mat'
        muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FCR', 'EPL', 'FDP2'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(1).Neut.Idx_25'-4500;
        
    case 4 %% D1 - Flx/Ext - Wrist Pronated 
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0032.mat'
        muscleName = {'APL', 'FPL', 'SUP', 'PTER', 'ECU', 'FDP2'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(1).Pro.Idx_25'-4500;
        
    case 5 %% D2 - Flx/Ext - Wrist Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0009.mat'
        muscleName = {'FDP2', 'FCU', 'SUP', 'PTER', 'ECU', 'FCR', 'ED2', 'EPL'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(2).Neut.Idx_25'-4500;

    case 6 %% D2 - Flx/Ext - Wrist Supinated
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0034.mat'
        muscleName = {'FDP2', 'FDS2', 'ED3', 'ED2','SUP', 'PTER', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(2).Sup.Idx_25'-4500;
        
    case 7 %% D4 - Flx/Ext - Wrist Neutral
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0011.mat'
        muscleName = {'FDS2', 'FDS4', 'ED3', 'EPL','SUP', 'PTER', 'ED4', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(4).Neut.Idx_25'-4500;
        
    case 8 %% D4 - Flx/Ext - Wrist Pronated 
        trialFile = 'RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0030.mat'
        muscleName = {'FDS2', 'FDS4', 'ED3', 'EPL','SUP', 'PTER', 'ED4', 'FCR'};
        cycleStarts = KinematicsIdx_Structured.Sbj(14).Digits(4).Pro.Idx_25'-4500;
        
    otherwise
        disp('Please enter a correct option or make a new action case!');
        return;
    end
    
    taskName = replace(trials({erase(trialFile, 'RNEL_')},:).sTrialType, '_', '-');
    
    % Identify Epochs based upon Case
    cycleStops = cycleStarts+Fs; %End Period for Data Epochs
    if length(cycleStops) ~= length(cycleStarts)
        cycleStarts = cycleStarts(1:length(cycleStops));
    end
    cycleLengths = (cycleStops - cycleStarts +1);

    for tIdx = 1:length(cycleLengths)
%     t_tmp = linspace(0,cycleLengths(tIdx)/Fs,cycleLengths(tIdx)); %
        t_tmp = linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
        epochTime{tIdx} = t_tmp; %linspace(cycleStarts(tIdx)/Fs,cycleStops(tIdx)/Fs,cycleLengths(tIdx));
        StartStops{actionCase}.Indices(tIdx,:) = [epochTime{tIdx}(1) epochTime{tIdx}(end)]
    end
    
    StartStops{actionCase}.trialname = taskName;
    
    
end

