%% UH3 EMG - Recruitment Curves - Quick Load Trial Stats %%
% The Purpose of this script is to quickly generate summary plots for the
% UH3 Aim 2: Reflexive EMG from stimulation Project.
% The Script utilizes the Stim Parameter .mat files (from Sensory Rig 1)
% and the trellis data.mat files (from Sensory Rig 2).  It will batch
% process 'n' number of files at once.  This should primarily be used for
% each step of the Recruitment Curve Process.  Please ensure that the
% target files are in the proper initial paths (RNELshare for stim params &
% the local datatank folder for trellis files).  Plots and Text Log will
% save in the designated Report Path for remote users to peruse.

Muscle_ID = {'Right GM', 'Right RF', 'Right VL', 'Right BF', 'Right SE', ...
    'Right TA', 'Right SO', 'Right LG', 'Left GM', 'Left RF', 'Left VL',...
    'Left BF', 'Left SE', 'Left VM', 'Left MG', 'Left LG'};


reportPath = 'R:\data_generated\human\uh3_stim\reports\testing';

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile('*nev','Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir(emgPathname, 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [erase(emgFilenames{f}, '.nev') '.mat'];
end

%% Initialize Analysis Parameters

% Filtering Settings:
fs = 30e3;
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 4;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  4;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs));

% Set Epoch Window for Stim-triggered Averaging
dtPre = 20e-3; % msec before stim onset 
dtPost = 80e-3; % msec after stim onset 


%% Begin Parsing each file:
for fnum = 1:length(emgFilenames)

    trialName(fnum,:) = erase(trialFilenames{fnum}, '.mat');
    %% Load Trial Stim Parameters
    if exist(fullfile(trialPathname,trialFilenames{fnum}), 'file')==0
        disp(['This Trial.mat DOES NOT EXIST - ' trialName(fnum,:)])
        fileExist{fnum}.trial = 0;
        fileExist{fnum}.name = trialName(fnum,:);
        continue;
    else
        disp(['Trial: ' erase(trialFilenames{fnum}, '.mat')]);
        fileExist{fnum}.trial = 1;
        fileExist{fnum}.name = trialName(fnum,:);
    end
       
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
    
    if isempty(trialInfo.Stim_params)
        disp(['This Trial.mat IS EMPTY - ' trialName(fnum,:)]);
        fileExist(fnum).trial = 2; %Where 2 equals exists but empty
        continue;
    end
    
    %% Identify Stim Events Channel
    spinalElec(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecs);
    nevStimCh(fnum) = trialInfo.Stim_params.NSChannels{1,1}(1);%find_stimevtchannel(spinalElec);
    
    % Set Pulse Width, Stim Duration, & Stim Frequency
    stimAmp(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecAmps);
    pulseWidth(fnum) = cell2mat(trialInfo.Stim_params.PulseWidth); %this is provided in ms
    stimDuration(fnum) = cell2mat(trialInfo.Stim_params.Duration)/1000; %this is provided in ms
    stimFrequency(fnum) = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
%     numStims(fnum) = stimDuration(fum)*stimFrequency(fnum);
    nSampsPre = floor(dtPre*fs);
    nSampsPost = floor(dtPost*fs);

    stimLength = (2*pulseWidth(fnum)/1000)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
    epochTimeVec = linspace(-dtPre, (stimLength/fs)+dtPost, nSampsPre+nSampsPost+stimLength);
    
    [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fnum}],nevStimCh(fnum));
    stims = floor(stimEvts{1}*fs);
    
  
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], 'raw' , 128+[1:32]);
    
    %Monopolar Recordings to Bipolar EMG
    for i = 1:2:size(analogData,1)
        epochRMS = analogData(i:i+1,:);
        emg(round(i/2),:) = diff(flipud(epochRMS));
    end
    
    emgF = filtfilt(b,a,emg');
    emgF = emgF';
    
    emgS = filtfilt(bS,aS,emgF')';
    % Epoch Windows 
    for i = 1:size(emg,1)
        asd = emgF(i,:);
        emgStim{i} = cell2mat(arrayfun(@(x) asd((stims(x)-(nSampsPre)):(stims(x)+(stimLength+nSampsPost-1))), 1:length(stims), 'UniformOutput',false)');
        qwe = emgS(i,:);
        emgStim_Sm{i} = cell2mat(arrayfun(@(x) asd((stims(x)-(nSampsPre)):(stims(x)+(stimLength+nSampsPost-1))), 1:length(stims), 'UniformOutput',false)');
         
    end
    
    %Plot Stim averages
    hF = figure;
    for c = 1:length(Muscle_ID)
        hFs(c) = subplot(2,round(length(Muscle_ID)/2),c);
        
     %Plot /Channel
     plot(epochTimeVec,emgStim{1,c});hold on;
     plot(epochTimeVec,mean(emgStim{1,c}),'k','LineWidth',1.2);
%      vline([0 stimLength/fs],'r-');
     title([Muscle_ID{c}]);
              
    end
    suplabel('Amplitude(uV)','y');
    suplabel('Time(sec)','x');
    
    Pix_SS = get(0,'screensize');
    hF.Position = Pix_SS;
    tmptit = suptitle([trialName(fnum,:) ' Stim Traces']);
    set(tmptit, 'Interpreter', 'none'); 
    linkaxes([hFs],'xy');
    hFs(1).XLim = [epochTimeVec(1) epochTimeVec(end)];

    disp('Making Trial Folder')
    mkdir(reportPath,trialName(fnum,:))
    savefig(hF,[fullfile(reportPath,trialName(fnum,:)) '\EpochedStimTraces']);
    saveas(hF,[fullfile(reportPath,trialName(fnum,:)) '\EpochedStimTraces.png']);

    %% Get RMS acros epochs
    epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
%     figure;cellfun(@plot,epochRMS);
    StimEnd = nSampsPre+stimLength;
    
    hF2 = figure;
    for c = 1:length(Muscle_ID)
              
     %Plot /Channel
     plot(epochTimeVec,epochRMS{1,c});hold on;
     %      vline([0 stimLength/fs],'r-');
     responseRMS(fnum,c) = mean(epochRMS{c}(StimEnd:end));
                  
    end
    xlabel('RMS Amplitude (uV)');
    ylabel('Time (Sec)');
    title([trialName(fnum,:) ' - Epoched RMS'],'Interpreter', 'none');
    vline(StimEnd/fs);
    legend(Muscle_ID);
    
    saveas(hF2,[fullfile(reportPath,trialName(fnum,:)) '\Epoched_RMS.png']);
   
    
    
    
    close all
    
    
    clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                dtPre dtPost fnum spinalElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency responseRMS trialName ...
                trialPathname trialFilenames
end

%% Plot Recruitment Curves

hF3 = figure;

subplot(2,3,1)
for m = 1:length(Muscle_ID)/2
    plot(stimAmp/1000,responseRMS(:,m)')
    hold on;
end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Amplitude (mA)')
    title('Recruitment Curve - Amplitude');
    
subplot(2,3,2)
for m = 1:length(Muscle_ID)/2
    scatter(pulseWidth,responseRMS(:,m)')
    hold on;
end
legend([Muscle_ID(1:length(Muscle_ID)/2)],'Location','northoutside', 'Orientation','horizontal');
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation PulseWidth (ms)')
    title('Recruitment Curve - PulseWidth');
    
subplot(2,3,3)
for m = 1:length(Muscle_ID)/2
    scatter(stimFrequency,responseRMS(:,m)')
    hold on;
end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Frequency (Hz)')
    title('Recruitment Curve - Frequency');
    
    subplot(2,3,4)
for m = (length(Muscle_ID)/2):length(Muscle_ID)
    plot(stimAmp/1000,responseRMS(:,m)')
    hold on;
end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Amplitude (mA)')
    title('Recruitment Curve - Amplitude');
    
subplot(2,3,5)
for m = (length(Muscle_ID)/2):length(Muscle_ID)
    scatter(pulseWidth,responseRMS(:,m)')
    hold on;
end
legend([Muscle_ID((length(Muscle_ID)/2)+1:length(Muscle_ID))],'Location','southoutside', 'Orientation','horizontal');
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation PulseWidth (ms)')
    title('Recruitment Curve - PulseWidth');
    
subplot(2,3,6)
for m = (length(Muscle_ID)/2):length(Muscle_ID)
    scatter(stimFrequency,responseRMS(:,m)')
    hold on;
end
    ylabel('(post-Stim) Mean RMS  (uV)')
    xlabel('Stimulation Frequency (Hz)')
    title('Recruitment Curve - Frequency');
    
Pix_SS = get(0,'screensize');
hF3.Position = Pix_SS;
    
saveas(hF3,[reportPath '\' erase(datestr(now),':') '_RecruitmentCurves.png']);
savefig(hF3,[reportPath '\' erase(datestr(now),':') '_RecruitmentCurves']);
    
    

varNames = {'TrialName','SpinalElectrodes', 'NEV_StimChannel', 'StimAmplitudes_mA', ...
    'PulseWidths_ms', 'StimDurations_s', 'StimFrequencies_Hz'};


TrialTable = table(string(trialName), spinalElec', nevStimCh', (stimAmp/1000)', pulseWidth',...
    stimDuration', stimFrequency','VariableNames', varNames)

 
writetable(TrialTable,[reportPath '\' erase(datestr(now),':') '_Summary.txt'],'Delimiter','tab')
writetable(TrialTable,[reportPath '\' erase(datestr(now),':') '_Summary.xls'])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

