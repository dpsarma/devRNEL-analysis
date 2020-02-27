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
for configCount = 1%1:numel(configFilenames)
    
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

    dtPre = 0.95; %seconds before stim onset to use | Due to Jitter, looking at 0.95 and not 1sec...
    dtPost = 1; %seconds after stim onset to use
    %% 
    % Set Epoch Window for Stim-triggered Averaging

    dtPre2 = 20/1000; %mseconds before stim onset to use
    dtPost2 = 50/1000; %mseconds after stim onset to use
    %% Start Analysis.  Need to Make Loop for ALL files.

    for fnum=1:length(emgFilenames) %1:2 were buggy |118 buggy |Afternoon, 122/180 buggy
%         tic
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
        if isempty(trialInfo.Stim_params)
            disp(['This Trial.mat IS EMPTY - ' trialName]);
            fileExist(fnum).trial = 2; %Where 2 equals exists but empty
            continue;
        end

        spinalElec = cell2mat(trialInfo.Stim_params.SpinalElecs);

        %Identify Stim Events Channel
        nevStimCh = trialInfo.Stim_params.NSChannels{1,1}(1);%find_stimevtchannel(spinalElec);

       
        %Load Data & Stim Events
        try
            [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fnum}],nevStimCh);
            if strcmp(trialName, 'LSP01_Ssn001_Set001_Blk001_Trl002')
                nevStimCh = 1;
            elseif strcmp(trialName, 'LSP01_Ssn006_Set023_Blk001_Trl036')
                nevStimCh = 53;
            end
        catch
            disp('An error occurred while retrieving stimEvts.');
            errFile{configCount, fnum} = emgFilenames{fnum};
            disp(['Recheck file: ' errFile(configCount, fnum)]);
            disp('Execution will continue.');
            continue;
        end
%         [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fnum}],nevStimCh);
        stims = stimEvts{1,1};
        
        [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], dtype , dchan); 

%         %% TESTING ONLY
%         timerCounter(fnum) = toc
%         continue; %% USE FOR TESTING "LOADING"   

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
            figure('visible', 'off');plot(tRef2,emgTmp);hold on;
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
            emgStimMean(c,:) = mean(emgPulseEpochs{1,c});
        end
    %% F. Plotting
    %%
        titString = sprintf('%s: %s - Channel %d, %dHz %dus %duA ', ssnName, trialName, spinalElec, stimFrequency, floor(stimPulse/fs*10^6), cell2mat(trialInfo.Stim_params.SpinalElecAmps));
    %% 
    % Plotting Full Stimulation Train "Response"


        % Stim Plot (16 channels, epoched)
        disp(['Plotting Full Stim Response - ' ssnString]); 
        figure('visible', 'off');
        for i = 1:size(emgStimFull,1)
            h(i) = subplot(size(emgStimFull,1)/2,2,i);
            plot(tRef,emgStimFull(i,:));
            vline(0,'r-', 'Stim Begin');
            vline(stimLength/fs,'r-', 'Stim Train End');
            axis([tRef(1) tRef(end) -inf inf]) %round(min(min(emgStimFull)),-3) round(max(max(emgStimFull)),-3)
            title(Muscle_ID(i));
            ylabel('EMG (uV)');
        end
        h(end-1).XLabel.String = 'Time (s)';
        h(end).XLabel.String = 'Time (s)';
        tmptit = suptitle(titString);
        set(tmptit, 'Interpreter', 'none'); 
        maximize;

        tmpDir = [ssnDir trialName '\'];
        mkdir(tmpDir)
        saveallopenfigs([tmpDir 'FullStimTrain-']);
    %% 
    % Plotting Stim-triggered "Average" Response
    %%
        % Single Stim-triggered EMG Plot (16 channels, epoched)
        disp(['Plotting Stim-triggered EMG Response - ' ssnString]); 
        figure('visible', 'off');
        for i = 1:size(emgStimRMS,1)
             h(i) = subplot(size(emgStimFull,1)/2,2,i);
            plot(tRef2,emgStimRMS(i,:));hold on;
            plot(tRef2,emgStimMean(i,:),'g');
            title(Muscle_ID(i));
            vline(0,'r-', 'Stim Pulse');
            vline(stimPulse/fs,'r-');
            axis([tRef2(1) tRef2(end) -inf inf])
            ylabel('EMG (uV)');
            legend('RMS (power)', 'mean','Location','NorthEast')
        end
        maximize;
        tmptit = suptitle({titString, 'StimTriggered RMS vs. Avg (green) EMG Responses'});
        set(tmptit, 'Interpreter', 'none');
        h(end-1).XLabel.String = 'Time (s)';
        h(end).XLabel.String = 'Time (s)';

        saveallopenfigs([ssnDir trialName '-']);


    %% 
    % Plotting Stim-triggered "Average" Response (Individually)
    %%
        % Stim Plots (All 16 channels, epoched)
        disp(['Plotting StimAvg by Muscle - ' ssnString]); 
        tmpDir = [ssnDir trialName '\stimEpochs\'];
        for i = 1:size(emgStimRMS,1)
            figure('visible', 'off');
            plot(tRef2,emgStimRMS(i,:));hold on;
            plot(tRef2,emgStimMean(i,:),'g');
            title(Muscle_ID(i));
            vline(0,'r-', 'Stim Pulse');
            vline(stimPulse/fs,'r-');
            axis([tRef2(1) tRef2(end) -inf inf])
            ylabel('EMG (uV)');
            xlabel('Time (s)');
            legend('RMS (power)', 'mean')
            title({titString, 'StimTriggered RMS vs. Avg (green) EMG Responses',Muscle_ID{i}}, 'Interpreter', 'none')
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
            figure('visible', 'off');
            hold on;
            plot(timeVec,emg(i,:));
            plot(timeVec,emgF(i,:),'g');
            title([titString Muscle_ID(i)], 'Interpreter', 'none');
            ylabel('EMG (uV)');
            xlabel('Time (s)');
            legend('raw','filt')
    %         ylim([round(min(min(emg)),-3) round(max(max(emg)),-3)])
            vline([stims]);
            maximize;
            saveallopenfigs([tmpDir Muscle_ID{i} '-']);
        end
    %     mkdir(['R:\data_generated\human\uh3_stim\LSP01\EMGresponses\scratch\' trialName '\raw\'])

    %%
%         timerCounter(fnum) = toc
    clearvars -except Muscle_ID Wp dtPost dtPre dtype emgPathname highCut ssnDir trialFilenames Norder a b dchan dtPost2 dtPre2 emgFilenames fs lowCut ssnName trialPathname fnum fileExist configCount configFilenames configPathname errFiles %timerCounter
    end

%     save([ssnName '-Timer.mat'],'timerCounter');
%     save([ssnName '-TrialExist.mat'],'fileExist');
    % save('UH3AnalTimer.mat','timerCounter');
    clearvars -except configCount configFilenames configPathname errFiles
end