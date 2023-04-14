%% Trial Summary UH3.  1/26/20
% This file goes through the trial files and summarizes what occured in
% each run.  Provide the data file and this reads and collates the headers.

addpath(genpath('R:\data_generated\human\uh3_stim\genSummary\code_tools')); 
path_datatank = 'R:\data_raw\human\uh3_stim\'; %'C:\DataTanks\2018\'; %'R:\data_raw\human\uh3_stim\';
path_stimparams = 'R:\data_raw\human\uh3_stim\'; %'R:\data_generated\human\uh3_stim\genSummary\'; %'R:\data_raw\human\uh3_stim\';

%% Experiment Session Parameters %'Limit Check f-2Hz PW-1ms' 'RC 500uA-5mA f-1Hz PW-200us'
% subjectName = 'LSP02b'; setName = 'Elec 1'; setDescrpt = 'Day 15 - 100Hz PW500us SD500ms - Sitting'; %file safe names only
% % subjectName = 'LSP05'; setName = 'All Days'; setDescrpt = 'All Trials with EMG'; %file safe names only
subjectName = 'LNP02'; setName = 'All Days'; setDescrpt = 'All Trials with EMG'; 
externalStimulator = 'no'; %no = NanoStims
chanPort = 128; % 128 or 256 based on which Grapevine Port (B or C)


% % reportPath = ['D:\FigRescources\UH3\' subjectName '\emgRecruitment_Summary\']; %['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];
reportPath = ['C:\data\LL_UH3\' subjectName '\emgRecruitment_Summary\'];
setPath = [reportPath setName '\'];
mkdir(setPath); %mkdir(reportPath);

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile([path_datatank subjectName '\*.nev'],'Pick files','MultiSelect', 'on'); %.nev or nf6
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir([path_stimparams subjectName '\'] , 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-10) '.mat']; %2014a friendly, no erase function (-4 for just file, -10 for _x# thing)
end



for fnum = 1:length(emgFilenames)
    
    disp(['Loading file ' num2str(fnum) ' of ' num2str(length(emgFilenames))]);
    trial(fnum).Name = erase(trialFilenames{fnum}, '.mat'); 
    
    %% Get Basic File Info
    try
        [ns_status, hFile] = ns_OpenFile([emgPathname emgFilenames{fnum}]); 
        [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
        ns_CloseFile(hFile);
    catch ME
        if (strcmp(ME.identifier,'Error using memmapfile/hChangeFilename'))
            disp(['Cannot Access File: ' trial(fnum).Name]);
            disp(['Not clear what the issue is with .cache'])
        else
            disp('Something SUPER weird happened w/.nev');
        end
        spinalElec{fnum} = nan;
        nevStimCh{fnum} = nan;
        stimAmp{fnum} = nan;
        pulseWidth(fnum) = nan; 
        stimDuration(fnum) = nan;
        stimFrequency(fnum) = nan;
        numStims(fnum) = nan;
        trialtype(fnum,:) = 'broken ';
        trial(fnum).Date = nan;
        trial(fnum).Time = nan';
        continue;
    end
    
    
    if isempty(hFile.Entity)
        disp([emgFilenames{fnum} ' is empty']);
        spinalElec{fnum} = nan;
        nevStimCh{fnum} = nan;
        stimAmp{fnum} = nan;
        pulseWidth(fnum) = nan; 
        stimDuration(fnum) = nan;
        stimFrequency(fnum) = nan;
        numStims(fnum) = nan;
        trialtype(fnum,:) = 'isempty';
        trial(fnum).Date = nan;
        trial(fnum).Time = nan';
        continue;
    end
    
    if contains(hFile.Entity(end).Label,'hifreq') %raw
        disp('File Contains EMG');
    else
        disp([emgFilenames{fnum} ' Does not have EMG']);
        spinalElec{fnum} = nan;
        nevStimCh{fnum} = nan;
        stimAmp{fnum} = nan;
        pulseWidth(fnum) = nan; 
        stimDuration(fnum) = nan;
        stimFrequency(fnum) = nan;
        numStims(fnum) = nan;
        trialtype(fnum,:) = 'no EMG ';
        trial(fnum).Date = nan;
        trial(fnum).Time = nan';
        continue;
    end
    
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
            disp('Something weird happened with trialInfo');
        end
        spinalElec{fnum} = nan;
        nevStimCh{fnum} = nan;
        stimAmp{fnum} = nan;
        pulseWidth(fnum) = nan; 
        stimDuration(fnum) = nan;
        stimFrequency(fnum) = nan;
        numStims(fnum) = nan;
        trialtype(fnum,:) = 'no .mat';
        trial(fnum).Date = nan;
        trial(fnum).Time = nan;
        continue;
    end
    
    if isempty(trialInfo.Stim_params)
        disp('Something weird happened with trialInfo');
        spinalElec{fnum} = nan;
        nevStimCh{fnum} = nan;
        stimAmp{fnum} = nan;
        pulseWidth(fnum) = nan; 
        stimDuration(fnum) = nan;
        stimFrequency(fnum) = nan;
        numStims(fnum) = nan;
        trialtype(fnum,:) = '0-trial';
        trial(fnum).Date = nan;
        trial(fnum).Time = nan;
        continue;
    end
    
    % Identify Stim Events Channel
        spinalElec{fnum} = (trialInfo.Stim_params(1).SpinalElecs);
        nevStimCh{fnum} = trialInfo.Stim_params(1).NSChannels{1,1}(1);   
        % Set Pulse Width, Stim Duration, & Stim Frequency
        stimAmp{fnum} = (trialInfo.Stim_params(1).SpinalElecAmps);
        pulseWidth(fnum) = cell2mat(trialInfo.Stim_params(1).PulseWidth); %provided in ms
        stimDuration(fnum) = cell2mat(trialInfo.Stim_params(1).Duration)/1000; %in s, provided in ms
        stimFrequency(fnum) = cell2mat(trialInfo.Stim_params(1).Frequency); %this is provided in Hz
        numStims(fnum) = stimDuration(fnum)*stimFrequency(fnum);
        
        switch trialInfo.Stim_params(1).trialType
            case 1
                trialtype(fnum,:) = 'sensory';
            otherwise
                trialtype(fnum,:) = 'other  ';
        end
        
         clearvars -except emgFilenames emgPathname reportPath  ...
                fnum spinalElec nevStimCh stimAmp pulseWidth trialtype...
                stimDuration stimFrequency trialName ...
                trialPathname trialFilenames trial subjectName setName ...
                setDescrpt setPath ...
                chanPort 
    close all
    
    disp('End File Processing');
end

%% Export Summary Table:
disp('Writing Summary Table');

varNames = {'TrialName','TrialType','SpinalElectrodes','Anode', 'NEV_StimChannel', 'StimAmplitudes_mA', 'AnodeAmp'...
    'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz', 'Record_Date', 'Record_Time'};

TrialNames = [{trial(:).Name}'];
TrialDates = [{trial(:).Date}'];
TrialTimes = [{trial(:).Time}'];

multpolar = (zeros(1,length(trial)));
multAmp = (zeros(1,length(trial)));
for i = 1:length(trial)
    try y(i) = cell2mat(spinalElec{i});
    catch
        disp('Hit a Nan')
        try y(i) = cell2mat(spinalElec(i));
        catch
            disp('Hit a Multiploar')
            y(i) = spinalElec{i}{:}(1);
            multpolar(i) = spinalElec{i}{:}(2);
        end
    end
    
    try yy(i) = cell2mat(stimAmp{i});
    catch
        disp('Hit a Nan')
        try yy(i) = cell2mat(stimAmp(i));
        catch
            disp('Hit a Multiploar')
            yy(i) = stimAmp{i}{:}(1);
            multAmp(i) = stimAmp{i}{:}(2);
        end
    end
end

TrialTable = table(TrialNames, trialtype, y',multpolar', nevStimCh', yy', multAmp', pulseWidth',...
    stimDuration', stimFrequency', TrialDates, TrialTimes, 'VariableNames', varNames)

 
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_' setDescrpt '_Summary.txt'],'Delimiter','tab');
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_' setDescrpt '_Summary.xls']);