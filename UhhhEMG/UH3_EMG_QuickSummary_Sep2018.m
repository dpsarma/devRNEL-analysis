%% UH3 EMG - Recruitment Curves - Quick Load Trial Stats - Testing Version %%
% The Purpose of this script is to quickly generate summary plots for the
% UH3 Aim 2: Reflexive EMG from stimulation Project.
% The Script utilizes the Stim Parameter .mat files (from Sensory Rig 1)
% and the trellis data.mat files (from Sensory Rig 2).  It will batch
% process 'n' number of files at once.  This should primarily be used for
% each step of the Recruitment Curve Process.  Please ensure that the
% target files are in the proper initial paths (RNELshare for stim params &
% the local datatank folder for trellis files).  Plots and Text Log will
% save in the designated Report Path for remote users to peruse.

%% FOR THE EXPERIMENTER:
% 1) Designate the Subject Name. 
% 2) Change Report Paths & Data Tank Path if needed. ('Stimparam'.mat files)
% 2) Ascribe Set (session/block, etc.).
% 3) Add a description (i.e. 1mA-6mA_Standing).
% 4) Set RMS Window and artifact settle buffer as needed.


addpath(genpath('R:\data_generated\human\uh3_stim\genSummary\code_tools')); 

%% Experiment Session Parameters
subjectName = 'DigitimerToo'; setName = 'Set_3'; setDescrpt = 'E1_1mA-12mA_Standing_tibial'; %file safe names only
% subjectName = 'LSP02'; setName = 'Set_1'; setDescrpt = 'E1_1mA-6mA_Standing'; %file safe names only
externalStimulator = 'yes'; %no = NanoStims
chanPort = 256; % 128 or 256 based on which Grapevine Port (B or C)
artifactBuffer = 0.015; % post stim artifact settle lag in s, so 0.01 is 10ms
rmsWindow = 0.05; % post stim RMS Window in s, so 0.05 is 50ms

% Muscles Used
Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right MG', 'Right SO lat', ...
    'Right TA', 'Right LG', 'Right SO'};

% Muscle_ID = {'Right GM', 'Right RF', 'Right VL', 'Right BF', 'Right SE', ...
%     'Right TA', 'Right SO', 'Right LG', 'Left GM', 'Left RF', 'Left VL',...
%     'Left BF', 'Left SE', 'Left VM', 'Left MG', 'Left LG'};

C = linspecer(length(Muscle_ID),'qualitative');

reportPath = ['R:\data_generated\human\uh3_stim\' subjectName '\emgRecruitment_Summary\'];
setPath = [reportPath setName '\'];
mkdir(setPath); %mkdir(reportPath);

%% Identify files to load
disp('Please Select the Ripple Data Folder');
[emgFilenames, emgPathname] = uigetfile(['C:\DataTanks\2018\' subjectName '\*.nev'],'Pick files','MultiSelect', 'on');
disp(['User selected ', fullfile(emgPathname, emgFilenames)]);

% disp('Please Select the Folder with the Stim Trial Info:');
trialPathname = uigetdir(['R:\data_generated\human\uh3_stim\genSummary\' subjectName '\'] , 'OpenLoop Stim Trial Info');

for f = 1:length(emgFilenames)
    trialFilenames{f} = [emgFilenames{f}(1:end-4) '.mat']; %2014a friendly, no erase function
end

%% Initialize Analysis Parameters

% Filtering Settings:
fs = 30e3;
lowCut = 75; %lowest frequency to pass
highCut = 7500; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

%Smoothing Settings:
Fsmooth = 10; %Lowpass cutoff frequency for smoothing the rectified signal 10 vs 50
NorderS =  2;
dt = 1/fs;
[bS,aS] = butter(NorderS,Fsmooth/(0.5*fs));

% Set Epoch Window for Stim-triggered Averaging
dtPre = 50e-3; % msec before stim onset 
dtPost = 100e-3; % msec after stim onset 
nSampsPre = floor(dtPre*fs);
nSampsPost = floor(dtPost*fs);

%% Begin Parsing each file:
for fnum = 1:length(emgFilenames)
    
    disp(['Loading file ' num2str(fnum) ' of ' num2str(length(emgFilenames))]);
    
    %% Get Basic File Info
    trial(fnum).Name = erase(trialFilenames{fnum}, '.mat'); 
    [ns_status, hFile] = ns_OpenFile([emgPathname emgFilenames{fnum}]); 
    [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
    
    trial(fnum).Date = join(string([nsFileInfo.Time_Month nsFileInfo.Time_Day nsFileInfo.Time_Year]),'-');
    trial(fnum).Time = join(string([nsFileInfo.Time_Hour nsFileInfo.Time_Min nsFileInfo.Time_Sec nsFileInfo.Time_MilliSec]),':');
    
    %% Load Data
    [analogData,timeVec] = read_continuousData([emgPathname emgFilenames{fnum}], 'raw' , chanPort+[1:length(Muscle_ID)*2]); %128 vs 256 | 8 vs. 16
    %% 'Data' -> Bipolar EMG
    for i = 1:2:size(analogData,1)
        bipolar = analogData(i:i+1,:);
        emg(round(i/2),:) = diff(flipud(bipolar));
    end
    %% Filtering (Denoising) & Smoothing
    emgF = filtfilt(b,a,emg')';  
    emgS = filtfilt(bS,aS,emgF')';
    
    disp('Filtering and Epoching....')
    
    %% Load Trial Stim Parameters
    trialInfo = load(fullfile(trialPathname,trialFilenames{fnum}));
        % Identify Stim Events Channel
        spinalElec(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecs);
        nevStimCh(fnum) = trialInfo.Stim_params.NSChannels{1,1}(1);   
        % Set Pulse Width, Stim Duration, & Stim Frequency
        stimAmp(fnum) = cell2mat(trialInfo.Stim_params.SpinalElecAmps);
        pulseWidth(fnum) = cell2mat(trialInfo.Stim_params.PulseWidth); %provided in ms
        stimDuration(fnum) = cell2mat(trialInfo.Stim_params.Duration)/1000; %in s, provided in ms
        stimFrequency(fnum) = cell2mat(trialInfo.Stim_params.Frequency); %this is provided in Hz
        numStims(fnum) = stimDuration(fnum)*stimFrequency(fnum);
    
    %% KAVEATS for Testing Data:
    
    if contains(trial(fnum).Name,'Set008')
        stimAmp(fnum) = 5900; %Correct for weird '76' fault
    end
    
    if contains(trial(fnum).Name,'Set012')
        stimAmp(fnum) = 5000 + cell2mat(trialInfo.Stim_params.SpinalElecAmps); %To go beyond 12mA
        if stimAmp(fnum) == 10200
            stimAmp(fnum) = 12000;
        elseif stimAmp(fnum) == 10400
            stimAmp(fnum) = 14000;
        elseif stimAmp(fnum) == 10600
            stimAmp(fnum) = 16000;
        elseif stimAmp(fnum) == 10800
            stimAmp(fnum) = 18000;
        elseif stimAmp(fnum) == 10800
            stimAmp(fnum) = 18000;
        elseif stimAmp(fnum) == 11000
            stimAmp(fnum) = 20000; 
        end
    end
    
    %% Load Stim Times
    switch externalStimulator
        case 'yes'
            [stimEvts] = read_digitalEvents([emgPathname emgFilenames{fnum}],1);
            stims = floor(cell2mat(stimEvts.timeStamp)*fs);
        otherwise
            [stimEvts,stchannels] = read_stimEvents([emgPathname emgFilenames{fnum}],nevStimCh(fnum));
            stims = floor(cell2mat(stimEvts.timeStamp)*fs);
    end
    
    stimLength = (pulseWidth(fnum)/1000)*fs + 2; %Pulse Width in ms, interstim interval is 2 ticks of clock (2/fs)
    trial(fnum).stimLength = stimLength;
    
    if length(stims) == 2*numStims(fnum)
        stimBegins = stims(1:2:end);
        stimEnds = stims(2:2:end);
    else
        stimBegins = stims;
    end
    %% Build Epoch Time Vector
    
    trial(fnum).epochTimeVec = linspace(-dtPre, (stimLength/fs) + dtPost, nSampsPre + nSampsPost + stimLength);
   
    %% Epoch data into Windows 
    for i = 1:size(emg,1)
        asd = emgF(i,:);
        emgStim{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
        qwe = emgS(i,:);
        emgStim_Sm{i} = cell2mat(arrayfun(@(x) asd((stimBegins(x)-nSampsPre):(stimBegins(x)+(stimLength+nSampsPost-1))), 1:length(stimBegins), 'UniformOutput',false)');
    end
    
    trial(fnum).emgStim = emgStim;
    trial(fnum).emgStim_Sm = emgStim_Sm;
    
    %% Get RMS acros epochs
    trial(fnum).epochRMS = cellfun(@rms,emgStim,'UniformOutput', false);
        
    %% Extract "Mean" Trace & Response
    stimEnd = nSampsPre + stimLength + artifactBuffer*fs;    
    for m = 1:length(Muscle_ID)
        % Mean Trace
        trial(fnum).meanTrace{m} = mean(emgStim{1,m});
        trial(fnum).meanTrace_Sm{m} = mean(emgStim_Sm{1,m});
        % Response
        trial(fnum).responseRMS(m,:) = mean(trial(fnum).epochRMS{m}(stimEnd:stimEnd+(rmsWindow*fs)));
        trial(fnum).responseP2P(m,:) = mean(trial(fnum).emgStim{m}(stimEnd:stimEnd+(rmsWindow*fs)));
    end

    %% Plot RMS Window for Verification
    disp('Plotting RMS Window for Verification...')
    
    hF2 = figure('visible', 'off');
        for c = 1:length(Muscle_ID)

         %%Plot /Channel
         plot(trial(fnum).epochTimeVec,trial(fnum).epochRMS{c},'Color',C(c,:));hold on;
              vline([0 stimLength/fs],'r-'); %vline(stimEnd/fs,'g-');

        end
        ylabel('RMS Amplitude (uV)');
        xlabel('Time (Sec)');
        title([trial(fnum).Name ' - Epoched RMS'],'Interpreter', 'none');
        vline([artifactBuffer artifactBuffer+rmsWindow], '-g');
        legend(Muscle_ID);    
    saveas(hF2,[setPath 'RMS_window-' trial(fnum).Name '.png']);
 
    %% clean variables for later plotting.
    
    clearvars -except emgFilenames emgPathname reportPath Muscle_ID fs b a bS aS ...
                dtPre dtPost fnum spinalElec nevStimCh stimAmp pulseWidth...
                stimDuration stimFrequency responseRMS trialName ...
                trialPathname trialFilenames trial subjectName setName ...
                setDescrpt setPath artifactBuffer rmsWindow ...
                externalStimulator nSampsPre nSampsPost chanPort C
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

 
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_Summary.txt'],'Delimiter','tab');
writetable(TrialTable,[reportPath '\' strrep(datestr(now),':','_') '_' setName '_Summary.xls']);

%% Plotting Stim Averages
for c = 1:length(Muscle_ID)
    disp('Plotting Stim Averages');
    
    hF(c) = figure; %('visible', 'off');
    tmpsup = suptitle([Muscle_ID{c} ' - ' setDescrpt ' -  Stim Traces']);
    set(tmpsup, 'Interpreter', 'none') 
    maximize;
%     set(hF(c),'Position', [1 -798 1080 1803]);
    subCount = 1;
    for fnum = 1:length(emgFilenames)
        hFs(fnum) = subplot(length(emgFilenames),1,fnum);%round(fnum/subcount)); evens by 2
%         hFs(fnum) = subplot(12,1,subCount);
        %Plot /Channel
        plot(trial(fnum).epochTimeVec,trial(fnum).emgStim{c});hold on;
        plot(trial(fnum).epochTimeVec,trial(fnum).meanTrace{c},'k','LineWidth',1.2);
        vline([0    trial(fnum).stimLength/fs],'r-');
        tmptit = title([trial(fnum).Name ' - ' num2str(stimAmp(fnum)/1000) 'mA, ' ...
            num2str(pulseWidth(fnum)) 'us, ' num2str(stimFrequency(fnum)) 'Hz, ' num2str(stimDuration(fnum)) 's']); 
        set(hFs(fnum),'FontSize',7);
        set(tmptit, 'Interpreter', 'none','FontSize',5);
        hFs(fnum).XLim = [trial(fnum).epochTimeVec(1) trial(fnum).epochTimeVec(end)];
        
        subCount = subCount + 1;
    end
    linkaxes([hFs],'xy');
    suplabel('Amplitude(uV)','y');
    suplabel('Time(sec)','x');
%     hF(c).Position = [317 49 568 1787];  
    
    disp('Saving to Set Folder')
    savefig(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName]);
    saveas(hF(c),[setPath Muscle_ID{c} '_EpochedStimTraces_' setName '.png']);
end

close all

%% Plotting Recruitment Curves
response = [trial(:).responseRMS]';
% response = [trial(:).responseP2P]';

switch length(Muscle_ID)
    case 16
        pRow = 2; pCol = 3;
        mCount = length(Muscle_ID)/2;
        pIdx = 2;
    case 8
        pRow = 1; pCol = 3;
        mCount = length(Muscle_ID);
        pIdx = 1;
    otherwise
        disp('How Would you like to Plot the Recruitment Curves?');
        return;
end

hF3 = figure; %maximize;
[sortedStim,I] = sort(stimAmp);
[sortedPulse,I2] = sort(pulseWidth);
[sortedFreq,I3] = sort(stimFrequency);

for pI = 1:pIdx
    
    subplot(pRow,pCol, 3*(pI-1)+1)
    for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
        yyA = smooth(sortedStim/1000,response(I,m),'sgolay');
        plot(sortedStim/1000,yyA,'Color',C(m,:),'LineWidth',3);
%         scatter(stimAmp/1000, response(I,m)');
        hold on;
    end
        ylabel('(post-Stim) Mean RMS  (uV)')
        xlabel('Stimulation Amplitude (mA)')
        title(['By Amplitude']);

    subplot(pRow,pCol, 3*(pI-1)+2)
    for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
%         scatter(pulseWidth,response(:,m)')
        yyP = smooth(sortedPulse,response(I2,m),'sgolay');
        plot(sortedPulse,yyP,'Color',C(m,:),'LineWidth',3);
        hold on;
    end
        ylabel('(post-Stim) Mean RMS  (uV)')
        xlabel('Stimulation PulseWidth (us)')
        title('By PulseWidth');
    %     
    subplot(pRow,pCol, 3*(pI-1)+3)
    for m = 1+mCount*(pI-1):mCount+mCount*(pI-1)
%         scatter(stimFrequency,response(:,m)')
        yyF = smooth(sortedFreq,response(I3,m),'sgolay');
        plot(sortedFreq,yyF,'Color',C(m,:),'LineWidth',3);
        hold on;
    end
        ylabel('(post-Stim) Mean RMS  (uV)')
        xlabel('Stimulation Frequency (Hz)')
        title('By Frequency');
        
    if pIdx == 1
        legend([Muscle_ID(1:mCount)],'Location','eastoutside', 'Orientation','vertical');
    elseif pIdx == 2
        legend([Muscle_ID((mCount)+1:length(Muscle_ID))],'Location','eastoutside', 'Orientation','vertical');
    end
%     
end
tmptit = suptitle([setDescrpt ' -  Recruitment Curves']);
set(tmptit, 'Interpreter', 'none');
Pix_SS = get(0,'screensize');
hF3.Position = Pix_SS;
    
saveas(hF3,[reportPath '\' strrep(datestr(now),':','_')  '_RecruitmentCurves' '_' setName '.png']);
savefig(hF3,[reportPath '\' strrep(datestr(now),':','_') '_RecruitmentCurves_' setName]);