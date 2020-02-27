%% UH3 EMG Response Report: Script A
% This live-script will attempt to generate a series of relevant plots to identify 
% EMG responses from each day's experiments (UH3-Aim2).
% 
% Need: Trellis Data Folder (30K recordings) & Stim Table Folder (trial info 
% from OpenStim) [Make Config File]
%% Load Config File & Initialize Filtering/Epoching Parameters
%%
disp('Please Select the Config File for the desired Session:')

load(uigetfile('*mat'));
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

%% 
% Set Filtering Parameters:

%%Filtering Parameters
if strcmp(dtype,'raw')
    fs = 30000;
elseif strcmp(dtype,'hi-res')
    fs = 2000;
end
    
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 4;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%% 
% Set Epoch Window for Full Stim Train

dtPre = 1; %seconds before stim onset to use
dtPost = 1; %seconds after stim onset to use
%% 
% Set Epoch Window for Stim-triggered Averaging

dtPre2 = 10/1000; %mseconds before stim onset to use
dtPost2 = 10/1000; %mseconds after stim onset to use
%% Start Analysis.  Need to Make Loop for ALL files.

for fnum=118:length(emgFilenames) %1:2 were buggy |118 buggy
    tic
%% A. Loading Data (in Meta Struct?)
%%
%    
    trialName = erase(trialFilenames{fnum}, '.mat');
    ssnString = [ssnName '_' trialName];
        
    %Load Trial Stim Parameters
    if exist(fullfile(trialPathname,trialFilenames{fnum}), 'file')==0
        disp(['This Trial.mat DOES NOT EXIST - ' trialName])
        fileExist{fnum}.trial = 0;
        fileExist{fnum}.name = trialName;
        continue;
    else
        disp(['Trial: ' erase(trialFilenames{fnum}, '.mat')]);
        fileExist{fnum}.trial = 1;
        fileExist{fnum}.name = trialName;
    end
        
        
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
    spinalElec = cell2mat(trialInfo.Stim_params.SpinalElecs);
    
    %Identify Stim Events Channel
    nevStimCh = find_stimevtchannel(spinalElec);
    
    if fnum ==2
        nevStimCh = 1; %This is ONLY for the silly Day-0 Files
    end
    
    %Load Data & Stim Events
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], dtype , dchan); 
    [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fnum}],nevStimCh);
    stims = stimEvts{1,1};

%% E. Differencing & Filtering
%%
    %Monopolar Recordings to Bipolar EMG
    emg = make_bipolar(analogData);
       
    %% Filtering 
    emgF = filtfilt(b,a,emg');
    emgF = emgF';
    
%% Stim Averaging
%%
    %Epoch Stim Data for Full Train
    [emgStimFull,tRef, stimLength] = epoch_stimEMG(emgF,stims, fs, [dtPre dtPost], 'full train'); 
    
    % Set Pulse Width, Stim Duration, & Stim Frequency
    pulseWidth = cell2mat(trialInfo.Stim_params.PulseWidth); %this is provided in ms
    stimDuration = cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
    stimFrequency = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
    numStims = stimDuration*stimFrequency;
    
    emgTmp = [];
    
    %Epoch Stim Data for Pulses (for Stim Averaging)
    tmpDir = [ssnDir trialName '\stimEpochs\'];
    mkdir(tmpDir)
%   mkdir(['R:\data_generated\human\uh3_stim\LSP01\EMGresponses\scratch\' trialName '\stimEpochs\'])
    for c = 1:size(emgF,1)
        display(['Extracting Epochs...' cell2mat(Muscle_ID(c))]);
        [emgTmp, tRef2, stimPulse] = epoch_stimEMG(emgF(c,:),stims,fs,[dtPre2 dtPost2],'pulses', pulseWidth, numStims);
        emgPulseEpochs{c} = emgTmp;
                
        %Plot /Channel
        figure;plot(tRef2,emgTmp);hold on;
        vline([0 stimPulse/fs],'r-');
        title([Muscle_ID{c} ': Stim Epochs']);

        saveallopenfigs([tmpDir Muscle_ID{c} '-']);
    end

  
    % RMS 'Average' (Power)
    for c=1:size(emgF,1)
        emgStimRMS(c,:) = rms(emgPulseEpochs{1,c});
    end
    
     % Mean 'Average'
    for c=1:size(emgF,1)
        emgStimMean(c,:) = rms(emgPulseEpochs{1,c});
    end
%% F. Plotting
%%
    titString = sprintf('Channel %d, %dHz %dus %dmA ', ssnName, trialName, spinalElec, stimFrequency, floor(stimPulse/fs*10^6), cell2mat(trialInfo.Stim_params.SpinalElecAmps));
%% 
% Plotting Full Stimulation Train "Response"

    
    % Stim Plot (16 channels, epoched)
    disp(['Plotting Full Stim Response - ' ssnString]); 
    figure;
    for i = 1:size(emgStimFull,1)
        h(i) = subplot(size(emgStimFull,1)/2,2,i);
        plot(tRef,emgStimFull(i,:));
        axis([tRef(1) tRef(end) round(min(min(emgStimFull)),-3) round(max(max(emgStimFull)),-3)])
        title(Muscle_ID(i));
        vline(0,'r-', 'Stim Pulse');
        vline(stimLength/fs,'r-');
        ylabel('EMG (uV)');
    end
    h(end-1).XLabel.String = 'Time (s)';
    h(end).XLabel.String = 'Time (s)';
    suptitle([ssnName ': ' titString])
    maximize;
    
    tmpDir = [ssnDir trialName '\'];
    mkdir(tmpDir)
    saveallopenfigs([tmpDir 'FullStimTrain-']);
%% 
% Plotting Stim-triggered "Average" Response
%%
    % Single Stim-triggered EMG Plot (16 channels, epoched)
    disp(['Plotting Stim-triggered EMG Response - ' ssnString]); 
    figure;
    for i = 1:size(emgStimRMS,1)
        subplot(8,2,i);
        plot(tRef2,emgStimRMS(i,:));hold on;
        plot(tRef2,emgStimMean(i,:),'g');
        title(Muscle_ID(i));
        vline(0,'r-', 'Stim Pulse');
        vline(stimPulse/fs,'r-');
        axis([tRef2(1) tRef2(end) -inf inf])
        ylabel('EMG (uV)');
        xlabel('Time (s)');
        legend('RMS (power)', 'mean','Location','EastOutside')
    end
    tmp = sprintf('Channel %d, %dHz %dus %dmA ', spinalElec, stimFrequency, floor(stimPulse/fs*10^6), cell2mat(trialInfo.Stim_params.SpinalElecAmps));
    suptitle({[ssnName ': ' tmp], 'StimTriggered RMS vs. Avg (green) EMG Responses'})
    maximize;
    
    saveallopenfigs([ssnDir trialName '-']);
    
    
%% 
% Plotting Stim-triggered "Average" Response (Individually)
%%
    % Stim Plots (All 16 channels, epoched)
    disp(['Plotting StimAvg by Muscle - ' ssnString]); 
    tmpDir = [ssnDir trialName '\stimEpochs\'];
    for i = 1:size(emgStimRMS,1)
        figure;
        plot(tRef2,emgStimRMS(i,:));hold on;
        plot(tRef2,emgStimMean(i,:),'g');
        title(Muscle_ID(i));
        vline(0,'r-', 'Stim Pulse');
        vline(stimPulse/fs,'r-');
        axis([tRef2(1) tRef2(end) -inf inf])
        ylabel('EMG (uV)');
        xlabel('Time (s)');
        legend('RMS (power)', 'mean')
        tmp = sprintf('Channel %d, %dHz %dus %dmA ', spinalElec, stimFrequency, floor(stimPulse/fs*10^6), cell2mat(trialInfo.Stim_params.SpinalElecAmps));
        title({[ssnName ': ' tmp], 'StimTriggered RMS vs. Avg (green) EMG Responses',Muscle_ID{i}})
        maximize;
        saveallopenfigs([tmpDir Muscle_ID{c} '-']);
    end
  
%% 
% Filter vs. Raw Plot (Full timeseries)
%%
    disp(['Plotting Filter vs. Raw Plot - ' ssnString]); 
    tmpDir = [ssnDir trialName '\raw\'];
    mkdir(tmpDir)    
    for i = 1:size(emg,1)
        figure;
        hold on;
        plot(timeVec,emg(i,:));
        plot(timeVec,emgF(i,:),'g');
        title([titString Muscle_ID(i)]);
        ylabel('EMG (uV)');
        xlabel('Time (s)');
        legend('raw','filt')
        ylim([round(min(min(emg)),-3) round(max(max(emg)),-3)])
        vline([stims(1) stims(end)]);
        saveallopenfigs([tmpDir Muscle_ID{i} '-']);
    end
%     mkdir(['R:\data_generated\human\uh3_stim\LSP01\EMGresponses\scratch\' trialName '\raw\'])
      
%%
    timerCounter(fnum) = toc
clearvars -except Muscle_ID Wp dtPost dtPre dtype emgPathname highCut ssnDir trialFilenames Norder a b dchan dtPost2 dtPre2 emgFilenames fs lowCut ssnName trialPathname fnum timerCounter  fileExist
end


save('UH3AnalTimer.mat','timerCounter');