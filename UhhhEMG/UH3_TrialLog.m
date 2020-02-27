%% Trial Summary UH3.  1/26/20
% This file goes through the trial files and summarizes what occured in
% each run.  Provide the data file and this reads and collates the headers.

addpath(genpath('R:\data_generated\human\uh3_stim\genSummary\code_tools')); 
path_datatank = 'R:\data_raw\human\uh3_stim\'; %'C:\DataTanks\2018\'; %'R:\data_raw\human\uh3_stim\';
path_stimparams = 'R:\data_raw\human\uh3_stim\'; %'R:\data_generated\human\uh3_stim\genSummary\'; %'R:\data_raw\human\uh3_stim\';

%% Experiment Session Parameters %'Limit Check f-2Hz PW-1ms' 'RC 500uA-5mA f-1Hz PW-200us'
% subjectName = 'LSP02b'; setName = 'Elec 1'; setDescrpt = 'Day 15 - 100Hz PW500us SD500ms - Sitting'; %file safe names only
subjectName = 'LSP05'; setName = 'Day12-All'; setDescrpt = 'RCs 1to8 plus PW-FQ mod on 16'; %file safe names only
externalStimulator = 'no'; %no = NanoStims
chanPort = 128; % 128 or 256 based on which Grapevine Port (B or C)


reportPath = ['D:\FigRescources\UH3\' subjectName '\emgRecruitment_Summary\']; %['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];
setPath = [reportPath setName '\'];
mkdir(setPath); %mkdir(reportPath);

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile([path_datatank subjectName '\*.nev'],'Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir([path_stimparams subjectName '\'] , 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-10) '.mat']; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end



for fnum = 1:length(emgFilenames)
    
    disp(['Loading file ' num2str(fnum) ' of ' num2str(length(emgFilenames))]);
       
    %% Get Basic File Info
    trial(fnum).Name = erase(trialFilenames{fnum}, '.mat'); 
    [ns_status, hFile] = ns_OpenFile([emgPathname emgFilenames{fnum}]); 
    
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    ns_CloseFile(hFile);
    
    trial(fnum).Date = join(string([nsFileInfo.Time_Month nsFileInfo.Time_Day nsFileInfo.Time_Year]),'-');
    trial(fnum).Time = join(string([nsFileInfo.Time_Hour nsFileInfo.Time_Min nsFileInfo.Time_Sec nsFileInfo.Time_MilliSec]),':');
    
    %% Load Trial Stim Parameters

    try
        trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
    catch ME
        if (strcmp(ME.identifier,'Error using load'))
            disp(['Unable to read file: ' trial(fnum).Name]);
            disp(['Most likely OpenLoop crashed and there is no .mat file'])
        else
            disp('Something weird happened');
        end
        spinalElec(fnum) = nan;
        nevStimCh(fnum) = nan;
        stimAmp(fnum) = nan;
        pulseWidth(fnum) = nan; 
        stimDuration(fnum) = nan;
        stimFrequency(fnum) = nan;
        numStims(fnum) = nan;
        continue;
    end
    
    
    
    % Identify Stim Events Channel
        spinalElec(fnum) = cell2mat(trialInfo.Stim_params(1).SpinalElecs);
        nevStimCh(fnum) = trialInfo.Stim_params(1).NSChannels{1,1}(1);   
        % Set Pulse Width, Stim Duration, & Stim Frequency
        stimAmp(fnum) = cell2mat(trialInfo.Stim_params(1).SpinalElecAmps);
        pulseWidth(fnum) = cell2mat(trialInfo.Stim_params(1).PulseWidth); %provided in ms
        stimDuration(fnum) = cell2mat(trialInfo.Stim_params(1).Duration)/1000; %in s, provided in ms
        stimFrequency(fnum) = cell2mat(trialInfo.Stim_params(1).Frequency); %this is provided in Hz
        numStims(fnum) = stimDuration(fnum)*stimFrequency(fnum);
        
        
         clearvars -except emgFilenames emgPathname reportPath  ...
                fnum spinalElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency trialName ...
                trialPathname trialFilenames trial subjectName setName ...
                setDescrpt setPath ...
                chanPort 
    close all
    
    disp('End File Processing');
end

%% Export Summary Table:
disp('Writing Summary Table');

varNames = {'TrialName','SpinalElectrodes', 'NEV_StimChannel', 'StimAmplitudes_mA', ...
    'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz', 'Record_Date', 'Record_Time'};

TrialNames = [{trial(:).Name}'];
TrialDates = [{trial(:).Date}'];
TrialTimes = [{trial(:).Time}'];

TrialTable = table(TrialNames, spinalElec', nevStimCh', (stimAmp/1000)', pulseWidth',...
    stimDuration', stimFrequency', TrialDates, TrialTimes, 'VariableNames', varNames)

 
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_' setDescrpt '_Summary.txt'],'Delimiter','tab');
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_' setDescrpt '_Summary.xls']);