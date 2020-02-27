%% UH3 EMG Response Report: Script A
% This live-script will attempt to generate a series of relevant plots to identify 
% EMG responses from each day's experiments (UH3-Aim2).
% 
% Need: Trellis Data Folder (30K recordings) & Stim Table Folder (trial info 
% from OpenStim) [Make Config File]
%% Load Config File & Initialize Filtering/Epoching Parameters
%%
disp('Please Select the Config Files for the desired Sessions:')

[configFilenames, configPathname] = uigetfile('*mat','Select Session Config files','MultiSelect', 'on');
errFile = {};
for configCount = 1%:numel(configFilenames)
    
    load(fullfile(configPathname,configFilenames{configCount}));
    % ExpType = {'Seated', 'Standing', 'Recruitment Curve'};

    disp(['Starting Analysis for: ' ssnName]);
    disp(['There are ' num2str(length(emgFilenames)) ' files from this session']);
    ssnDir = ['D:\Figures\UH3\LSPO1\EMGResponses\' ssnName '\'];
    mkdir(ssnDir); 

    % Set Data File Type (nf3/hi-res or ns5/raw)
    dtype = 'raw';
    % dtype = 'hi-res';
    disp(['Analyzing ' dtype ' data from Channels ' num2str(dchan(1)) ' to ' num2str(dchan(end))]);

    %% Should Add Environment Variables for Saving.... Future ADD!!!!


    %% Start Analysis.  Need to Make Loop for ALL files.

    for fnum=1:length(emgFilenames) %1:2 were buggy |118 buggy |Afternoon, 122/180 buggy
    %% A. Loading Data (in Meta Struct?)
    %%
    %    
        tmp = string(erase(trialFilenames{fnum}, '.mat'));
        trialName(fnum) = tmp;
        ssnString = [ssnName '_' trialName];

        %Load Trial Stim Parameters
        if exist(fullfile(trialPathname,trialFilenames{fnum}), 'file')==0
            disp(['This Trial.mat DOES NOT EXIST - ' trialName]);
            fileExist{fnum}.trial = 0;
            fileExist{fnum}.name = trialName;
            continue;
        else
            disp(['Trial: ' erase(trialFilenames{fnum}, '.mat')]);
            fileExist{fnum}.trial = 1;
            fileExist{fnum}.name = trialName;
        end


        trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
        if isempty(trialInfo.Stim_params)
            disp(['This Trial.mat IS EMPTY - ' trialName]);
%             fileExist(fnum).trial = 2; %Where 2 equals exists but empty
            continue;
        end

        spinalElec(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecs);
        stimAmplitude(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecAmps);
        stimPulseWidth(fnum) = cell2mat(trialInfo.Stim_params.PulseWidth);
        stimDuration(fnum) = cell2mat(trialInfo.Stim_params.Duration);
        stimFrequency(fnum) = cell2mat(trialInfo.Stim_params.Frequency);
        
        %Identify Stim Events Channel
        nevStimCh(fnum) = trialInfo.Stim_params.NSChannels{1,1}(1);%find_stimevtchannel(spinalElec);

        %Load Data & Stim Events
    
    disp([ssnString]);

    clearvars -except trialName spinalElec nevStimCh stimAmplitude stimPulseWidth stimDuration stimFrequency Muscle_ID emgPathname ssnDir trialFilenames emgFilenames ssnName trialPathname fnum fileExist configCount configFilenames configPathname errFiles
    end

    trialTable = table(trialName',spinalElec',nevStimCh',stimAmplitude',stimPulseWidth',stimDuration',stimFrequency');
    save([ssnName '-trialTable.mat'],'trialTable');
    writetable(trialTable, [ssnName '-trialTable.txt']);
    clearvars -except configCount configFilenames configPathname errFiles 
end
clear all
close all