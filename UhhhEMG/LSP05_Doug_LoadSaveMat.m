%% Loading Script for LSP05 Data (for Doug)

path_datatank = 'R:\data_raw\human\uh3_stim\'; %'C:\DataTanks\2018\'; %'R:\data_raw\human\uh3_stim\';
path_stimparams = 'R:\data_raw\human\uh3_stim\'; %'R:\data_generated\human\uh3_stim\genSummary\'; %'R:\data_raw\human\uh3_stim\';
path_datagen = 'R:\data_generated\human\uh3_stim\LSP05\';

%% Experiment Session Parameters %'Limit Check f-2Hz PW-1ms' 'RC 500uA-5mA f-1Hz PW-200us'
% subjectName = 'LSP02b'; setName = 'Elec 1'; setDescrpt = 'Day 15 - 100Hz PW500us SD500ms - Sitting'; %file safe names only
subjectName = 'LSP05'; setName = 'Day8-Set_1'; setDescrpt = 'E9-16RC_1-6mA_500us_1Hz'; %file safe names only

chanPort = 128; % 128 or 256 based on which Grapevine Port (B or C)

% Muscles Used in Channel Order
Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
    'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF', 'Left VL',...
    'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile([path_datatank subjectName '\*.nev'],'Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir([path_stimparams subjectName '\'] , 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-10) '.mat']; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end

%% Begin Parsing each file:
for fnum = 1:length(emgFilenames)
    
    disp(['Loading file ' num2str(fnum) ' of ' num2str(length(emgFilenames))]);
    
    %% Get Basic File Info
    trial(fnum).Name = erase(trialFilenames{fnum}, '.mat'); 
    [ns_status, hFile] = ns_OpenFile([emgPathname emgFilenames{fnum}]); 
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    ns_CloseFile(hFile);
    
    trial(fnum).Date = join(string([nsFileInfo.Time_Month nsFileInfo.Time_Day nsFileInfo.Time_Year]),'-');
    trial(fnum).Time = join(string([nsFileInfo.Time_Hour nsFileInfo.Time_Min nsFileInfo.Time_Sec nsFileInfo.Time_MilliSec]),':');
    
    %% Load Data
    [trial(fnum).Data,trial(fnum).timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], 'raw' , chanPort+[1:length(Muscle_ID)*2]); %128 vs 256 | 8 vs. 16

    %% Load Trial Stim Parameters
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
        % Identify Stim Events Channel
        
        trial(fnum).nevStimCh = trialInfo.Stim_params(1).NSChannels{1,1}(1);   
        % Set Pulse Width, Stim Duration, & Stim Frequency
        trial(fnum).pulseWidth = cell2mat(trialInfo.Stim_params(1).PulseWidth); %provided in ms
        trial(fnum).stimDuration = cell2mat(trialInfo.Stim_params(1).Duration)/1000; %in s, provided in ms
        trial(fnum).stimFrequency = cell2mat(trialInfo.Stim_params(1).Frequency); %this is provided in Hz
        trial(fnum).numStims = trial(fnum).stimDuration*trial(fnum).stimFrequency;
end

save([path_datagen setName '_' setDescrpt '.mat'],'trial', '-v7.3');